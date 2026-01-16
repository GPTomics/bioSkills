#!/bin/bash
# Screen FASTQ files for contamination
# Usage: ./screen_samples.sh <fastq_dir> <output_dir>

INPUT_DIR="${1:-.}"
OUTPUT_DIR="${2:-screen_results}"
THREADS="${3:-8}"

echo "=== Contamination Screening ==="
echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"

mkdir -p "$OUTPUT_DIR"

FASTQ_FILES=$(find "$INPUT_DIR" -name "*.fastq.gz" -o -name "*.fq.gz" | head -20)
N_FILES=$(echo "$FASTQ_FILES" | wc -l)

if [[ -z "$FASTQ_FILES" ]]; then
    echo "No FASTQ files found in $INPUT_DIR"
    exit 1
fi

echo "Found $N_FILES FASTQ files to screen"

for fq in $FASTQ_FILES; do
    sample=$(basename "$fq" | sed 's/.fastq.gz//' | sed 's/.fq.gz//')
    echo -e "\n=== Screening: $sample ==="

    fastq_screen \
        --threads "$THREADS" \
        --outdir "$OUTPUT_DIR" \
        "$fq"
done

echo -e "\n=== Aggregating Results ==="
if command -v multiqc &> /dev/null; then
    multiqc "$OUTPUT_DIR" -o "${OUTPUT_DIR}/multiqc" -f
    echo "MultiQC report: ${OUTPUT_DIR}/multiqc/multiqc_report.html"
else
    echo "MultiQC not installed, skipping aggregation"
fi

echo -e "\n=== Summary ==="
echo "Results directory: $OUTPUT_DIR"
echo "Individual reports: ${OUTPUT_DIR}/*_screen.html"

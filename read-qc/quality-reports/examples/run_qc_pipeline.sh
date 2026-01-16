#!/bin/bash
# Quality control pipeline for FASTQ files
# Usage: ./run_qc_pipeline.sh <input_dir> <output_dir>

INPUT_DIR="${1:-.}"
OUTPUT_DIR="${2:-qc_results}"
THREADS="${3:-4}"

echo "=== Quality Control Pipeline ==="
echo "Input: $INPUT_DIR"
echo "Output: $OUTPUT_DIR"
echo "Threads: $THREADS"

mkdir -p "$OUTPUT_DIR/fastqc"
mkdir -p "$OUTPUT_DIR/multiqc"

echo -e "\n=== Running FastQC ==="
fastqc -t "$THREADS" -o "$OUTPUT_DIR/fastqc" "$INPUT_DIR"/*.fastq.gz

if [[ $? -ne 0 ]]; then
    echo "FastQC failed"
    exit 1
fi

echo -e "\n=== Running MultiQC ==="
multiqc "$OUTPUT_DIR/fastqc" -o "$OUTPUT_DIR/multiqc" -f

if [[ $? -ne 0 ]]; then
    echo "MultiQC failed"
    exit 1
fi

echo -e "\n=== Summary ==="
echo "FastQC reports: $OUTPUT_DIR/fastqc/"
echo "MultiQC report: $OUTPUT_DIR/multiqc/multiqc_report.html"

N_SAMPLES=$(ls "$OUTPUT_DIR/fastqc"/*_fastqc.html 2>/dev/null | wc -l)
echo "Samples processed: $N_SAMPLES"

echo -e "\n=== Done ==="

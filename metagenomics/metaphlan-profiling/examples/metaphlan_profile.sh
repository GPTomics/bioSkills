#!/bin/bash
# Reference: MetaPhlAn 4.1+, Bowtie2 2.5.3+ | Verify API if version differs
# Profile multiple samples on ONE pinned database index, cache each mapping for cheap
# re-profiling, then merge. MetaPhlAn percentages are cell fractions (taxonomic abundance),
# never to be merged with Kraken/Bracken read fractions.
set -euo pipefail

READS_DIR="fastq"
OUTPUT_DIR="metaphlan_output"
NPROC=8
INDEX="mpa_vJun23_CHOCOPhlAnSGB_202403"   # pin the DB version; it is a batch variable across a study

mkdir -p "$OUTPUT_DIR/profiles" "$OUTPUT_DIR/mapout"

for fq in "${READS_DIR}"/*.fastq.gz; do
    sample=$(basename "$fq" .fastq.gz)
    echo "Processing ${sample}..."
    # Pre-4.2 builds use --bowtie2out instead of --mapout; check `metaphlan --help`.
    metaphlan "$fq" \
        --input_type fastq \
        --index "$INDEX" \
        --nproc "$NPROC" \
        --output_file "${OUTPUT_DIR}/profiles/${sample}_profile.txt" \
        --mapout "${OUTPUT_DIR}/mapout/${sample}.map.bz2"
done

# All inputs must share the same --index or rows mismatch silently.
echo "Merging profiles..."
merge_metaphlan_tables.py "${OUTPUT_DIR}"/profiles/*_profile.txt > "${OUTPUT_DIR}/merged_abundance.txt"

# Species rows only (exclude the t__ SGB tier). These are cell fractions, not read fractions.
echo "Top species across all samples:"
grep "s__" "${OUTPUT_DIR}/merged_abundance.txt" | grep -v "t__" | head -20
echo "Results saved to ${OUTPUT_DIR}/merged_abundance.txt"

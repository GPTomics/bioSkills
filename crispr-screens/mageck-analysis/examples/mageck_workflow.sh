#!/bin/bash
# Complete MAGeCK CRISPR screen analysis

set -e

# Configuration
LIBRARY="library.csv"
OUTPUT_PREFIX="screen_analysis"

# Input FASTQ files
DAY0="Day0.fastq.gz"
CTRL1="Control_rep1.fastq.gz"
CTRL2="Control_rep2.fastq.gz"
TREAT1="Treatment_rep1.fastq.gz"
TREAT2="Treatment_rep2.fastq.gz"

echo "=== Step 1: Count sgRNAs ==="
mageck count \
    -l "$LIBRARY" \
    -n "$OUTPUT_PREFIX" \
    --sample-label Day0,Control1,Control2,Treatment1,Treatment2 \
    --fastq "$DAY0" "$CTRL1" "$CTRL2" "$TREAT1" "$TREAT2" \
    --norm-method median

echo "=== Step 2: MAGeCK Test ==="
mageck test \
    -k "${OUTPUT_PREFIX}.count.txt" \
    -t Treatment1,Treatment2 \
    -c Control1,Control2 \
    -n "${OUTPUT_PREFIX}_test" \
    --norm-method median

echo "=== Results Summary ==="
echo "Gene summary: ${OUTPUT_PREFIX}_test.gene_summary.txt"
echo "sgRNA summary: ${OUTPUT_PREFIX}_test.sgrna_summary.txt"

# Count significant hits
NEG_HITS=$(awk -F'\t' 'NR>1 && $8<0.05' "${OUTPUT_PREFIX}_test.gene_summary.txt" | wc -l)
POS_HITS=$(awk -F'\t' 'NR>1 && $14<0.05' "${OUTPUT_PREFIX}_test.gene_summary.txt" | wc -l)

echo "Negative selection hits (FDR<0.05): $NEG_HITS"
echo "Positive selection hits (FDR<0.05): $POS_HITS"

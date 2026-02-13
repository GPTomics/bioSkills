#!/bin/bash
# Reference: GATK 4.5+, bedtools 2.31+ | Verify API if version differs
# CNV calling with CNVkit

set -euo pipefail

TUMOR_BAM=${1:-tumor.bam}
NORMAL_BAM=${2:-normal.bam}
TARGETS=${3:-targets.bed}
REFERENCE=${4:-reference.fa}
OUTPUT_DIR=${5:-cnvkit_output}

mkdir -p $OUTPUT_DIR

echo "=== CNVkit Analysis ==="

# Build reference from normal
echo "Building reference..."
cnvkit.py batch $TUMOR_BAM \
    --normal $NORMAL_BAM \
    --targets $TARGETS \
    --fasta $REFERENCE \
    --output-reference ${OUTPUT_DIR}/reference.cnn \
    --output-dir $OUTPUT_DIR

# Call CNVs
echo "Calling CNVs..."
cnvkit.py call ${OUTPUT_DIR}/*.cns \
    -o ${OUTPUT_DIR}/calls.cns

# Generate plots
echo "Generating plots..."
cnvkit.py scatter ${OUTPUT_DIR}/*.cnr \
    -s ${OUTPUT_DIR}/*.cns \
    -o ${OUTPUT_DIR}/scatter.pdf

cnvkit.py diagram ${OUTPUT_DIR}/*.cnr \
    -s ${OUTPUT_DIR}/*.cns \
    -o ${OUTPUT_DIR}/diagram.pdf

echo "Running quality metrics..."
SAMPLE=$(basename $TUMOR_BAM .bam)
cnvkit.py metrics ${OUTPUT_DIR}/${SAMPLE}.cnr -s ${OUTPUT_DIR}/${SAMPLE}.cns

echo "Gene-level report..."
# 0.2: CNVkit default gain/loss log2 threshold (2^0.2 ~ 15% copy number change)
cnvkit.py genemetrics ${OUTPUT_DIR}/${SAMPLE}.cnr -s ${OUTPUT_DIR}/${SAMPLE}.cns \
    --threshold 0.2 -o ${OUTPUT_DIR}/${SAMPLE}.genemetrics.tsv

echo "Results in: $OUTPUT_DIR"

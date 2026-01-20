#!/bin/bash
# TF footprinting with TOBIAS

set -euo pipefail

BAM=${1:-atac.bam}
PEAKS=${2:-peaks.bed}
GENOME=${3:-genome.fa}
MOTIFS=${4:-motifs.jaspar}
OUTPUT_DIR=${5:-tobias_output}

mkdir -p $OUTPUT_DIR

echo "=== TOBIAS Footprinting ==="

# Step 1: Correct Tn5 bias
echo "Correcting Tn5 bias..."
TOBIAS ATACorrect \
    --bam $BAM \
    --genome $GENOME \
    --peaks $PEAKS \
    --outdir $OUTPUT_DIR \
    --cores 8

# Step 2: Calculate footprint scores
echo "Calculating footprint scores..."
TOBIAS FootprintScores \
    --signal ${OUTPUT_DIR}/*_corrected.bw \
    --regions $PEAKS \
    --output ${OUTPUT_DIR}/footprints.bw \
    --cores 8

# Step 3: Bind detection
echo "Detecting TF binding..."
TOBIAS BINDetect \
    --motifs $MOTIFS \
    --signals ${OUTPUT_DIR}/footprints.bw \
    --genome $GENOME \
    --peaks $PEAKS \
    --outdir ${OUTPUT_DIR}/bindetect \
    --cores 8

echo "Results in: ${OUTPUT_DIR}/bindetect/"

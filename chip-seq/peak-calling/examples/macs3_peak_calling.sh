#!/bin/bash

CHIP_BAM="chip.sorted.bam"
INPUT_BAM="input.sorted.bam"
OUTPUT_DIR="peaks"
SAMPLE_NAME="experiment"
GENOME="hs"

mkdir -p $OUTPUT_DIR

macs3 callpeak \
    -t $CHIP_BAM \
    -c $INPUT_BAM \
    -f BAM \
    -g $GENOME \
    -n ${SAMPLE_NAME}_narrow \
    --outdir $OUTPUT_DIR \
    -q 0.05 \
    -B --SPMR  # q=0.05 is MACS default FDR; use 0.01 for stricter, 0.1 for exploratory.

macs3 callpeak \
    -t $CHIP_BAM \
    -c $INPUT_BAM \
    -f BAM \
    -g $GENOME \
    -n ${SAMPLE_NAME}_broad \
    --outdir $OUTPUT_DIR \
    --broad \
    --broad-cutoff 0.1 \
    -B --SPMR  # broad-cutoff=0.1 is MACS default for linking subpeaks; use 0.05 for stricter boundaries.

echo "Narrow peaks: $(wc -l < ${OUTPUT_DIR}/${SAMPLE_NAME}_narrow_peaks.narrowPeak)"
echo "Broad peaks: $(wc -l < ${OUTPUT_DIR}/${SAMPLE_NAME}_broad_peaks.broadPeak)"

sort -k8,8nr ${OUTPUT_DIR}/${SAMPLE_NAME}_narrow_peaks.narrowPeak | head -10

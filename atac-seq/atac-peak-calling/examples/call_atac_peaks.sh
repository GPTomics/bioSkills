#!/bin/bash
# ATAC-seq peak calling with MACS3

set -euo pipefail

BAM=${1:-atac.bam}
OUTPUT_PREFIX=${2:-atac_peaks}
GENOME_SIZE=${3:-hs}

echo "=== ATAC-seq Peak Calling ==="

# ATAC-specific parameters:
# --nomodel: Don't build shifting model
# --shift -100 --extsize 200: Center on Tn5 cut sites
# -q 0.01: FDR threshold

macs3 callpeak \
    -t $BAM \
    -f BAMPE \
    -g $GENOME_SIZE \
    -n $OUTPUT_PREFIX \
    --nomodel \
    --shift -100 \
    --extsize 200 \
    -q 0.01 \
    --keep-dup all \
    -B \
    --SPMR

echo "Peaks: ${OUTPUT_PREFIX}_peaks.narrowPeak"
echo "Peak count: $(wc -l < ${OUTPUT_PREFIX}_peaks.narrowPeak)"

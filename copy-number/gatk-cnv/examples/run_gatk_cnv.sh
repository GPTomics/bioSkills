#!/bin/bash
# GATK CNV calling

set -euo pipefail

BAM=${1:-sample.bam}
PON=${2:-cnv_pon.hdf5}
INTERVALS=${3:-targets.interval_list}
OUTPUT_PREFIX=${4:-gatk_cnv}

echo "=== GATK CNV Analysis ==="

# Collect read counts
gatk CollectReadCounts \
    -I $BAM \
    -L $INTERVALS \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O ${OUTPUT_PREFIX}.counts.hdf5

# Denoise with PON
gatk DenoiseReadCounts \
    -I ${OUTPUT_PREFIX}.counts.hdf5 \
    --count-panel-of-normals $PON \
    --standardized-copy-ratios ${OUTPUT_PREFIX}.standardizedCR.tsv \
    --denoised-copy-ratios ${OUTPUT_PREFIX}.denoisedCR.tsv

# Model segments
gatk ModelSegments \
    --denoised-copy-ratios ${OUTPUT_PREFIX}.denoisedCR.tsv \
    --output-prefix $OUTPUT_PREFIX \
    -O .

# Call copy ratios
gatk CallCopyRatioSegments \
    -I ${OUTPUT_PREFIX}.cr.seg \
    -O ${OUTPUT_PREFIX}.called.seg

# Plot
gatk PlotDenoisedCopyRatios \
    --standardized-copy-ratios ${OUTPUT_PREFIX}.standardizedCR.tsv \
    --denoised-copy-ratios ${OUTPUT_PREFIX}.denoisedCR.tsv \
    --sequence-dictionary reference.dict \
    --output-prefix $OUTPUT_PREFIX \
    -O plots

echo "Results: ${OUTPUT_PREFIX}.called.seg"

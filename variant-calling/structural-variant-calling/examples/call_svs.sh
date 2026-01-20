#!/bin/bash
# Structural variant calling with Manta and Delly

set -euo pipefail

BAM=${1:-sample.bam}
REFERENCE=${2:-reference.fa}
OUTPUT_DIR=${3:-sv_results}

mkdir -p $OUTPUT_DIR

echo "=== Structural Variant Calling ==="

# Manta (better for insertions/deletions)
echo "Running Manta..."
configManta.py \
    --bam $BAM \
    --referenceFasta $REFERENCE \
    --runDir ${OUTPUT_DIR}/manta

${OUTPUT_DIR}/manta/runWorkflow.py -m local -j 8

# Delly (better for inversions, translocations)
echo "Running Delly..."
delly call \
    -g $REFERENCE \
    -o ${OUTPUT_DIR}/delly_raw.bcf \
    $BAM

delly filter \
    -f germline \
    -o ${OUTPUT_DIR}/delly_filtered.bcf \
    ${OUTPUT_DIR}/delly_raw.bcf

bcftools view ${OUTPUT_DIR}/delly_filtered.bcf -Oz -o ${OUTPUT_DIR}/delly.vcf.gz

# Merge calls
echo "Merging SV calls..."
SURVIVOR merge <(echo "${OUTPUT_DIR}/manta/results/variants/diploidSV.vcf.gz"; echo "${OUTPUT_DIR}/delly.vcf.gz") \
    1000 2 1 1 0 50 \
    ${OUTPUT_DIR}/merged_svs.vcf

echo "Results:"
echo "  Manta: ${OUTPUT_DIR}/manta/results/variants/diploidSV.vcf.gz"
echo "  Delly: ${OUTPUT_DIR}/delly.vcf.gz"
echo "  Merged: ${OUTPUT_DIR}/merged_svs.vcf"

#!/bin/bash
# Run DeepVariant on a sample BAM file

set -euo pipefail

BAM=${1:-sample.bam}
REFERENCE=${2:-reference.fa}
OUTPUT_PREFIX=${3:-deepvariant_output}
MODEL_TYPE=${4:-WGS}
THREADS=${5:-8}

echo "=== DeepVariant Variant Calling ==="
echo "BAM: $BAM"
echo "Reference: $REFERENCE"
echo "Model: $MODEL_TYPE"

# Run DeepVariant
docker run -v "${PWD}:/data" google/deepvariant:1.6.0 \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=${MODEL_TYPE} \
    --ref=/data/${REFERENCE} \
    --reads=/data/${BAM} \
    --output_vcf=/data/${OUTPUT_PREFIX}.vcf.gz \
    --output_gvcf=/data/${OUTPUT_PREFIX}.g.vcf.gz \
    --num_shards=${THREADS}

# Index outputs
bcftools index -t ${OUTPUT_PREFIX}.vcf.gz
bcftools index -t ${OUTPUT_PREFIX}.g.vcf.gz

# Generate statistics
bcftools stats ${OUTPUT_PREFIX}.vcf.gz > ${OUTPUT_PREFIX}_stats.txt

echo "=== Complete ==="
echo "VCF: ${OUTPUT_PREFIX}.vcf.gz"
echo "gVCF: ${OUTPUT_PREFIX}.g.vcf.gz"
echo "Stats: ${OUTPUT_PREFIX}_stats.txt"

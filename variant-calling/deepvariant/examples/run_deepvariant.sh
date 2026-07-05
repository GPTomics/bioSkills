#!/bin/bash
# Reference: DeepVariant 1.6.1+, bcftools 1.19+ | Verify API if version differs
# Run DeepVariant on a sample BAM. Input BAM must be sorted, indexed, and
# duplicate-marked -- do NOT run BQSR first (it lowers DeepVariant accuracy).

set -euo pipefail

BAM=${1:-sample.bam}
REFERENCE=${2:-reference.fa}
OUTPUT_PREFIX=${3:-deepvariant_output}
MODEL_TYPE=${4:-WGS}          # WGS|WES|PACBIO|ONT_R104|HYBRID_PACBIO_ILLUMINA; must match the instrument
THREADS=${5:-8}               # shards for the CPU-bound make_examples stage

echo "=== DeepVariant: ${MODEL_TYPE} mode ==="
echo "BAM: $BAM  Reference: $REFERENCE"

docker run -v "${PWD}:/data" google/deepvariant:1.6.1 \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=${MODEL_TYPE} \
    --ref=/data/${REFERENCE} \
    --reads=/data/${BAM} \
    --output_vcf=/data/${OUTPUT_PREFIX}.vcf.gz \
    --output_gvcf=/data/${OUTPUT_PREFIX}.g.vcf.gz \
    --num_shards=${THREADS}

bcftools index -t ${OUTPUT_PREFIX}.vcf.gz
bcftools index -t ${OUTPUT_PREFIX}.g.vcf.gz

# DeepVariant output is already CNN-filtered (FILTER = PASS/RefCall); do NOT
# apply GATK hard filters or VQSR. Stats only for QC (Ti/Tv ~2.0-2.1 for WGS).
bcftools stats ${OUTPUT_PREFIX}.vcf.gz > ${OUTPUT_PREFIX}_stats.txt

echo "=== Complete ==="
echo "VCF: ${OUTPUT_PREFIX}.vcf.gz"
echo "gVCF: ${OUTPUT_PREFIX}.g.vcf.gz"
echo "Stats: ${OUTPUT_PREFIX}_stats.txt"

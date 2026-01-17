#!/bin/bash
# Comprehensive variant annotation pipeline
set -euo pipefail

INPUT=$1
OUTPUT_PREFIX=$2
VEP_CACHE=${3:-$HOME/.vep}
SNPEFF_DB=${4:-GRCh38.105}

echo "=== Variant Annotation Pipeline ==="
echo "Input: $INPUT"
echo "Output prefix: $OUTPUT_PREFIX"

# Option 1: VEP annotation
echo "Running VEP..."
vep -i $INPUT \
    -o ${OUTPUT_PREFIX}_vep.vcf \
    --vcf \
    --cache --offline --dir_cache $VEP_CACHE \
    --assembly GRCh38 \
    --everything \
    --pick \
    --fork 4 \
    --force_overwrite

bgzip -f ${OUTPUT_PREFIX}_vep.vcf
tabix -p vcf ${OUTPUT_PREFIX}_vep.vcf.gz

# Option 2: SnpEff annotation
echo "Running SnpEff..."
snpEff ann \
    -v \
    -stats ${OUTPUT_PREFIX}_snpeff_stats.html \
    $SNPEFF_DB \
    $INPUT > ${OUTPUT_PREFIX}_snpeff.vcf

bgzip -f ${OUTPUT_PREFIX}_snpeff.vcf
tabix -p vcf ${OUTPUT_PREFIX}_snpeff.vcf.gz

# Filter high impact variants (VEP output)
echo "Filtering high impact variants..."
bcftools view -i 'INFO/CSQ~"HIGH"' \
    ${OUTPUT_PREFIX}_vep.vcf.gz \
    -Oz -o ${OUTPUT_PREFIX}_high_impact.vcf.gz
tabix -p vcf ${OUTPUT_PREFIX}_high_impact.vcf.gz

# Summary
echo "=== Summary ==="
echo "Total variants: $(bcftools view -H $INPUT | wc -l)"
echo "High impact: $(bcftools view -H ${OUTPUT_PREFIX}_high_impact.vcf.gz | wc -l)"
echo ""
echo "Output files:"
echo "  ${OUTPUT_PREFIX}_vep.vcf.gz - VEP annotated"
echo "  ${OUTPUT_PREFIX}_snpeff.vcf.gz - SnpEff annotated"
echo "  ${OUTPUT_PREFIX}_high_impact.vcf.gz - High impact only"

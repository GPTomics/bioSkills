#!/bin/bash
set -euo pipefail

INPUT_VCF=$1
OUTPUT_PREFIX=${2:-filtered}
REFERENCE=${3:-reference.fa}

echo "=== Splitting SNPs and Indels ==="
bcftools view -v snps "$INPUT_VCF" -Oz -o "${OUTPUT_PREFIX}_snps_raw.vcf.gz"
bcftools view -v indels "$INPUT_VCF" -Oz -o "${OUTPUT_PREFIX}_indels_raw.vcf.gz"

echo "=== Filtering SNPs ==="
bcftools filter -i '
    QUAL >= 30 &&
    INFO/DP >= 10 &&
    (INFO/QD >= 2.0 || INFO/QD = ".") &&
    (INFO/FS <= 60.0 || INFO/FS = ".") &&
    (INFO/MQ >= 40.0 || INFO/MQ = ".") &&
    (INFO/MQRankSum >= -12.5 || INFO/MQRankSum = ".") &&
    (INFO/ReadPosRankSum >= -8.0 || INFO/ReadPosRankSum = ".") &&
    (INFO/SOR <= 3.0 || INFO/SOR = ".")
' "${OUTPUT_PREFIX}_snps_raw.vcf.gz" -Oz -o "${OUTPUT_PREFIX}_snps_filtered.vcf.gz"

echo "=== Filtering Indels ==="
bcftools filter -i '
    QUAL >= 30 &&
    INFO/DP >= 10 &&
    (INFO/QD >= 2.0 || INFO/QD = ".") &&
    (INFO/FS <= 200.0 || INFO/FS = ".") &&
    (INFO/ReadPosRankSum >= -20.0 || INFO/ReadPosRankSum = ".") &&
    (INFO/SOR <= 10.0 || INFO/SOR = ".")
' "${OUTPUT_PREFIX}_indels_raw.vcf.gz" -Oz -o "${OUTPUT_PREFIX}_indels_filtered.vcf.gz"

echo "=== Merging filtered variants ==="
bcftools concat "${OUTPUT_PREFIX}_snps_filtered.vcf.gz" "${OUTPUT_PREFIX}_indels_filtered.vcf.gz" | \
    bcftools sort -Oz -o "${OUTPUT_PREFIX}_all_filtered.vcf.gz"
bcftools index -t "${OUTPUT_PREFIX}_all_filtered.vcf.gz"

echo "=== Statistics ==="
echo "Before filtering:"
bcftools stats "$INPUT_VCF" | grep -E "^SN|^TSTV" | head -20

echo ""
echo "After filtering:"
bcftools stats "${OUTPUT_PREFIX}_all_filtered.vcf.gz" | grep -E "^SN|^TSTV" | head -20

TSTV=$(bcftools stats "${OUTPUT_PREFIX}_all_filtered.vcf.gz" | grep "^TSTV" | cut -f5)
echo ""
echo "Ti/Tv ratio: $TSTV (expected ~2.1 for WGS, ~2.8-3.3 for exomes)"

echo ""
echo "Output: ${OUTPUT_PREFIX}_all_filtered.vcf.gz"

#!/bin/bash
# Genotype imputation with Beagle

set -euo pipefail

VCF=${1:-input.vcf.gz}
REF_PANEL=${2:-reference.vcf.gz}
MAP=${3:-plink.map}
OUTPUT=${4:-imputed}

echo "=== Genotype Imputation with Beagle ==="

java -Xmx16g -jar beagle.jar \
    gt=$VCF \
    ref=$REF_PANEL \
    map=$MAP \
    out=$OUTPUT \
    nthreads=8

# Index output
bcftools index -t ${OUTPUT}.vcf.gz

# Filter by imputation quality
bcftools view -i 'INFO/DR2>0.8' ${OUTPUT}.vcf.gz -Oz -o ${OUTPUT}_filtered.vcf.gz

echo "Imputed: ${OUTPUT}.vcf.gz"
echo "Filtered (DR2>0.8): ${OUTPUT}_filtered.vcf.gz"

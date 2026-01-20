#!/bin/bash
# Haplotype phasing with SHAPEIT5

set -euo pipefail

VCF=${1:-input.vcf.gz}
MAP=${2:-genetic_map.txt}
OUTPUT=${3:-phased.vcf.gz}

echo "=== Haplotype Phasing with SHAPEIT5 ==="

# Phase common variants
shapeit5_phase_common \
    --input $VCF \
    --map $MAP \
    --output ${OUTPUT%.vcf.gz}_common.bcf \
    --thread 8

# Phase rare variants (optional, requires phased common)
shapeit5_phase_rare \
    --input $VCF \
    --scaffold ${OUTPUT%.vcf.gz}_common.bcf \
    --map $MAP \
    --output $OUTPUT \
    --thread 8

echo "Phased output: $OUTPUT"

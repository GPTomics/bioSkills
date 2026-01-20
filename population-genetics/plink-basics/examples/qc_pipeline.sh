#!/bin/bash
# PLINK QC pipeline for population genetics
# Usage: ./qc_pipeline.sh <input_vcf> <output_prefix>

INPUT="${1}"
PREFIX="${2:-qc_output}"

if [[ -z "$INPUT" ]]; then
    echo "Usage: $0 <input.vcf.gz> [output_prefix]"
    exit 1
fi

echo "=== PLINK QC Pipeline ==="
echo "Input: $INPUT"
echo "Output prefix: $PREFIX"

echo -e "\n=== Step 1: Convert VCF to PLINK ==="
plink2 --vcf "$INPUT" --double-id --make-bed --out "${PREFIX}_raw"

echo -e "\n=== Step 2: Initial Statistics ==="
echo "Samples: $(wc -l < ${PREFIX}_raw.fam)"
echo "Variants: $(wc -l < ${PREFIX}_raw.bim)"

plink2 --bfile "${PREFIX}_raw" --missing --out "${PREFIX}_raw"

echo -e "\n=== Step 3: Apply QC Filters ==="
# Thresholds follow PLINK best practices for GWAS; adjust for rare variant studies
plink2 --bfile "${PREFIX}_raw" \
    --maf 0.01 \      # Minor allele freq; typical 0.01-0.05 depending on sample size
    --geno 0.05 \     # Max 5% missing genotypes per variant
    --mind 0.05 \     # Max 5% missing genotypes per sample
    --hwe 1e-6 \      # Hardy-Weinberg p-value; stringent to catch genotyping errors
    --make-bed --out "$PREFIX"

echo -e "\n=== Step 4: Final Statistics ==="
echo "Samples after QC: $(wc -l < ${PREFIX}.fam)"
echo "Variants after QC: $(wc -l < ${PREFIX}.bim)"

INITIAL_SNPS=$(wc -l < "${PREFIX}_raw.bim")
FINAL_SNPS=$(wc -l < "${PREFIX}.bim")
echo "Variant retention: $(echo "scale=1; $FINAL_SNPS * 100 / $INITIAL_SNPS" | bc)%"

INITIAL_SAMP=$(wc -l < "${PREFIX}_raw.fam")
FINAL_SAMP=$(wc -l < "${PREFIX}.fam")
echo "Sample retention: $(echo "scale=1; $FINAL_SAMP * 100 / $INITIAL_SAMP" | bc)%"

echo -e "\n=== Output Files ==="
ls -lh "${PREFIX}".{bed,bim,fam}

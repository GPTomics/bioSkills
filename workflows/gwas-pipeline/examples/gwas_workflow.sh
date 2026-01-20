#!/bin/bash
# Complete GWAS workflow with PLINK2
set -e

INPUT_VCF="genotypes.vcf.gz"
PHENO_FILE="phenotypes.txt"
OUTDIR="gwas_results"

mkdir -p ${OUTDIR}

echo "=== GWAS Pipeline ==="

# === Step 1: Import and Initial QC ===
echo "=== Step 1: Data Import and QC ==="

# Convert VCF to PLINK
plink2 --vcf ${INPUT_VCF} \
    --make-bed \
    --out ${OUTDIR}/raw

echo "Initial variants: $(wc -l < ${OUTDIR}/raw.bim)"
echo "Initial samples: $(wc -l < ${OUTDIR}/raw.fam)"

# Sample QC: remove high missing rate
plink2 --bfile ${OUTDIR}/raw \
    --mind 0.05 \
    --make-bed \
    --out ${OUTDIR}/sample_qc

# Variant QC
plink2 --bfile ${OUTDIR}/sample_qc \
    --geno 0.05 \
    --maf 0.01 \
    --hwe 1e-6 \
    --make-bed \
    --out ${OUTDIR}/qc

echo "After QC variants: $(wc -l < ${OUTDIR}/qc.bim)"
echo "After QC samples: $(wc -l < ${OUTDIR}/qc.fam)"

# === Step 2: LD Pruning ===
echo "=== Step 2: LD Pruning ==="

plink2 --bfile ${OUTDIR}/qc \
    --indep-pairwise 50 5 0.2 \
    --out ${OUTDIR}/ld_prune

echo "Independent variants: $(wc -l < ${OUTDIR}/ld_prune.prune.in)"

plink2 --bfile ${OUTDIR}/qc \
    --extract ${OUTDIR}/ld_prune.prune.in \
    --make-bed \
    --out ${OUTDIR}/pruned

# === Step 3: Population Structure ===
echo "=== Step 3: PCA ==="

plink2 --bfile ${OUTDIR}/pruned \
    --pca 10 \
    --out ${OUTDIR}/pca

# === Step 4: Association Testing ===
echo "=== Step 4: Association Testing ==="

# With PCA covariates (columns 3-12 are PC1-PC10)
plink2 --bfile ${OUTDIR}/qc \
    --pheno ${PHENO_FILE} \
    --covar ${OUTDIR}/pca.eigenvec \
    --covar-col-nums 3-12 \
    --glm hide-covar \
    --out ${OUTDIR}/gwas

# === Step 5: Extract Results ===
echo "=== Step 5: Processing Results ==="

# Find result file
result_file=$(ls ${OUTDIR}/gwas.*.glm.* 2>/dev/null | head -1)

if [ -f "$result_file" ]; then
    # Genome-wide significant
    awk 'NR==1 || $12 < 5e-8' "$result_file" > ${OUTDIR}/significant_5e8.txt
    sig_count=$(tail -n +2 ${OUTDIR}/significant_5e8.txt | wc -l)

    # Suggestive
    awk 'NR==1 || $12 < 1e-5' "$result_file" > ${OUTDIR}/suggestive_1e5.txt
    sug_count=$(tail -n +2 ${OUTDIR}/suggestive_1e5.txt | wc -l)

    echo ""
    echo "=== Results Summary ==="
    echo "Genome-wide significant (p < 5e-8): ${sig_count}"
    echo "Suggestive (p < 1e-5): ${sug_count}"
else
    echo "Warning: No result file found"
fi

echo ""
echo "=== GWAS Complete ==="
echo "Results directory: ${OUTDIR}/"
echo "  - QC'd data: ${OUTDIR}/qc.{bed,bim,fam}"
echo "  - PCA: ${OUTDIR}/pca.eigenvec"
echo "  - Association: ${OUTDIR}/gwas.*.glm.*"

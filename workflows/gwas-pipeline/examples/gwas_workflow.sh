#!/bin/bash
# Reference: PLINK 2.0 alpha 5+ | Verify API if version differs
# Complete GWAS workflow with PLINK2
set -e

INPUT_VCF="genotypes.vcf.gz"
PHENO_FILE="phenotypes.txt"
OUTDIR="gwas_results"

mkdir -p ${OUTDIR}

echo "=== GWAS Pipeline ==="

# === Step 1: Import and Initial QC ===
echo "=== Step 1: Data Import and QC ==="

# Convert VCF to PLINK (load --pheno now so PHENO1 is in the .fam for the controls-only HWE gate)
plink2 --vcf ${INPUT_VCF} \
    --pheno ${PHENO_FILE} \
    --make-bed \
    --out ${OUTDIR}/raw

echo "Initial variants: $(wc -l < ${OUTDIR}/raw.bim)"
echo "Initial samples: $(wc -l < ${OUTDIR}/raw.fam)"

# Variant missingness FIRST (own invocation), so samples are not dropped for soon-removed variants.
plink2 --bfile ${OUTDIR}/raw --geno 0.05 --make-bed --out ${OUTDIR}/var_qc

# Sample missingness, then KING-robust relatedness (structure-robust, IBD-free).
plink2 --bfile ${OUTDIR}/var_qc --mind 0.05 --king-cutoff 0.0884 --make-bed --out ${OUTDIR}/samp_qc

# HWE in CONTROLS ONLY (mid-p): a true non-additive risk variant deviates in cases, so a
# case-inclusive HWE test would remove real signal. Then MAF.
plink2 --bfile ${OUTDIR}/samp_qc --keep-if "PHENO1 == control" --hwe 1e-6 midp \
    --write-snplist --out ${OUTDIR}/hwe_pass
plink2 --bfile ${OUTDIR}/samp_qc --maf 0.01 --extract ${OUTDIR}/hwe_pass.snplist \
    --make-bed --out ${OUTDIR}/qc

echo "After QC variants: $(wc -l < ${OUTDIR}/qc.bim)"
echo "After QC samples: $(wc -l < ${OUTDIR}/qc.fam)"

# === Step 2: LD Pruning ===
echo "=== Step 2: LD Pruning ==="

# Exclude long-range-LD/inversion regions (MHC, 8p23.1, 17q21.31, LCT) FIRST -- their internal r2 is
# high and real, so a window prune keeps them and a PC then tracks the inversion, not ancestry.
plink2 --bfile ${OUTDIR}/qc --exclude range longrange_ld.txt --make-bed --out ${OUTDIR}/noLR
plink2 --bfile ${OUTDIR}/noLR --indep-pairwise 50 5 0.1 --out ${OUTDIR}/ld_prune

echo "Independent variants: $(wc -l < ${OUTDIR}/ld_prune.prune.in)"

plink2 --bfile ${OUTDIR}/noLR \
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
    --glm firth-fallback hide-covar \
    --out ${OUTDIR}/gwas

# === Step 5: Extract Results ===
echo "=== Step 5: Processing Results ==="

# Find result file
result_file=$(ls ${OUTDIR}/gwas.*.glm.* 2>/dev/null | head -1)

if [ -f "$result_file" ]; then
    # Select the P column by HEADER, not a fixed index (firth-fallback shifts columns).
    awk 'NR==1{for(i=1;i<=NF;i++) if($i=="P") p=i; print; next} $p<5e-8' \
        "$result_file" > ${OUTDIR}/significant_5e8.txt
    sig_count=$(tail -n +2 ${OUTDIR}/significant_5e8.txt | wc -l)

    awk 'NR==1{for(i=1;i<=NF;i++) if($i=="P") p=i; print; next} $p<1e-5' \
        "$result_file" > ${OUTDIR}/suggestive_1e5.txt
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

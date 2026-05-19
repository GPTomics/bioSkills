#!/bin/bash
# Reference: PRSice-2 2.3.5+, plink2 2.00a3+, ldsc 1.0+ | Verify CLI flags if version differs
#
# PRSice-2 clumping + thresholding workflow. C+T is the LEGACY baseline -- prefer
# LDpred2-auto (Prive 2020), SBayesRC (Zheng 2024), or PRS-CSx for multi-ancestry.
# PRSice-2 is retained for: published clinical scores defined via C+T (e.g., PRS313 / Mavaddat 2019);
# quick exploratory analyses; comparison baseline against Bayesian methods.

set -euo pipefail

GWAS="${GWAS:-gwas_summary.txt}"
TARGET="${TARGET:-target_genotypes}"
PHENO="${PHENO:-phenotype.txt}"
COV="${COV:-covariates.txt}"
OUT="${OUT:-prs_output}"

if [ ! -f "$GWAS" ]; then
    echo "Required: GWAS sumstats with columns SNP CHR BP A1 A2 BETA SE P"
    exit 1
fi
if [ ! -f "${TARGET}.bed" ]; then
    echo "Required: target plink files (.bed/.bim/.fam)"
    exit 1
fi

# Pre-flight: bivariate LDSC for sample-overlap detection (REQUIRED before any PRS evaluation)
echo "=== Pre-flight: sample overlap check ==="
echo "Run bivariate LDSC manually if a discovery sumstats and target sumstats both exist."
echo "Threshold: |gcov_int| > 0.05 with target n >= 1000 indicates substantial overlap."
echo "If confirmed, do NOT proceed with evaluation -- rebuild disjoint splits."
echo

# Strand-ambiguous SNPs: PRSice-2 default drops them (a/t c/g pairs).
# Clumping: r^2 < 0.1, 250 kb window -- conventional.
# p-value thresholds: standard bar-level sweep.

echo "=== PRSice-2 C+T (legacy baseline) ==="
PRSice_linux \
    --base "$GWAS" \
    --target "$TARGET" \
    --snp SNP --chr CHR --bp BP --A1 A1 --A2 A2 --pvalue P --beta BETA \
    --clump-kb 250 --clump-r2 0.1 \
    --bar-levels 5e-8,1e-5,1e-3,0.01,0.05,0.1,0.5,1 \
    --fastscore --all-score \
    --thread 4 \
    --out "$OUT"

echo "PRSice-2 complete. Output:"
echo "  ${OUT}.summary -- best threshold + R^2"
echo "  ${OUT}.all_score -- scores at each threshold"
echo "  ${OUT}.prsice -- per-SNP info"

# Validated run with phenotype + covariates (PCs in evaluation regression)
# CRITICAL: PCs must be from TEST cohort, not discovery cohort
if [ -f "$PHENO" ] && [ -f "$COV" ]; then
    echo
    echo "=== PRSice-2 with phenotype + covariates (PCs in eval) ==="
    PRSice_linux \
        --base "$GWAS" \
        --target "$TARGET" \
        --pheno "$PHENO" \
        --cov "$COV" \
        --cov-col @PC[1-10],Age,Sex \
        --binary-target T \
        --snp SNP --chr CHR --bp BP --A1 A1 --A2 A2 --pvalue P --beta BETA \
        --clump-kb 250 --clump-r2 0.1 \
        --thread 4 \
        --out "${OUT}_validated"
fi

echo
echo "=== Reviewer-grade follow-up ==="
echo "1. Switch to LDpred2-auto or SBayesRC for production use (Bayesian shrinkage)."
echo "2. Run EraSOR / bivariate LDSC for sample-overlap detection."
echo "3. Exclude HLA region (chr6 28-34 Mb) from main PRS; model classical HLA separately."
echo "4. Apply ancestry-conditional Z normalization using TEST cohort PCs (not discovery)."
echo "5. Transform PRS percentile to absolute risk via external incidence curve."
echo "6. Report Hingorani 2023 caveats: HR per SD = 1.3-1.7 is comparable to family history."
echo "7. PRS-RS 22-item checklist: priority items 13 (cohort independence), 16 (confounders),"
echo "   19 (absolute risk), 21 (ancestry composition of validation)."

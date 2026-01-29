#!/bin/bash
# PRSice-2 workflow for polygenic risk score calculation
#
# Prerequisites:
#   - PRSice_linux binary in PATH
#   - GWAS summary statistics
#   - Target genotypes in plink format (bed/bim/fam)
#
# Example data:
#   - Download GWAS from GWAS Catalog: https://www.ebi.ac.uk/gwas/
#   - Download 1000 Genomes for testing: https://www.internationalgenome.org/

# Input files (customize for your data)
GWAS="gwas_summary.txt"
TARGET="target_genotypes"  # plink prefix (bed/bim/fam)
PHENO="phenotype.txt"      # Optional: FID IID pheno
COV="covariates.txt"       # Optional: FID IID PC1 PC2 ... Age Sex
OUT="prs_output"

# Check inputs exist
if [ ! -f "$GWAS" ]; then
    echo "GWAS file not found: $GWAS"
    echo "Create example GWAS file:"
    echo "SNP CHR BP A1 A2 BETA SE P"
    echo "rs12345 1 10000 A G 0.05 0.01 1e-8"
    exit 1
fi

if [ ! -f "${TARGET}.bed" ]; then
    echo "Target genotypes not found: ${TARGET}.bed"
    echo "Provide plink format files (bed/bim/fam)"
    exit 1
fi

# Basic PRSice-2 run
# --clump-kb 250: clump SNPs within 250kb
# --clump-r2 0.1: remove SNPs with r2 > 0.1
# --bar-levels: p-value thresholds to test
# --fastscore: output scores at bar-levels only (faster)
# --all-score: output scores for all samples

echo "Running PRSice-2..."
PRSice_linux \
    --base "$GWAS" \
    --target "$TARGET" \
    --snp SNP \
    --chr CHR \
    --bp BP \
    --A1 A1 \
    --A2 A2 \
    --pvalue P \
    --beta BETA \
    --clump-kb 250 \
    --clump-r2 0.1 \
    --bar-levels 5e-8,1e-5,1e-3,0.01,0.05,0.1,0.5,1 \
    --fastscore \
    --all-score \
    --thread 4 \
    --out "$OUT"

echo "PRSice-2 complete. Output files:"
echo "  ${OUT}.summary - Best threshold and R2"
echo "  ${OUT}.all_score - Scores at all thresholds"
echo "  ${OUT}.prsice - Per-SNP and threshold info"

# With phenotype and covariates for validation
if [ -f "$PHENO" ] && [ -f "$COV" ]; then
    echo ""
    echo "Running with phenotype and covariates..."
    PRSice_linux \
        --base "$GWAS" \
        --target "$TARGET" \
        --pheno "$PHENO" \
        --cov "$COV" \
        --cov-col @PC[1-10],Age,Sex \
        --binary-target T \
        --snp SNP \
        --chr CHR \
        --bp BP \
        --A1 A1 \
        --A2 A2 \
        --pvalue P \
        --beta BETA \
        --clump-kb 250 \
        --clump-r2 0.1 \
        --thread 4 \
        --out "${OUT}_validated"

    echo "Validated output: ${OUT}_validated.summary"
fi

echo ""
echo "Next steps:"
echo "1. Check ${OUT}.summary for best p-value threshold"
echo "2. Use ${OUT}.all_score for downstream analysis"
echo "3. Normalize scores and compute percentiles"

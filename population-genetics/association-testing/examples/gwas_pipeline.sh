#!/bin/bash
# GWAS pipeline with PCA correction
# Usage: ./gwas_pipeline.sh <plink_prefix> <pheno_file> <output_prefix>

BFILE="${1}"
PHENO="${2}"
PREFIX="${3:-gwas_results}"

if [[ -z "$BFILE" ]] || [[ -z "$PHENO" ]]; then
    echo "Usage: $0 <plink_prefix> <phenotype_file> [output_prefix]"
    exit 1
fi

echo "=== GWAS Pipeline ==="
echo "Input data: $BFILE"
echo "Phenotypes: $PHENO"
echo "Output prefix: $PREFIX"

echo -e "\n=== Step 1: Calculate PCs ==="
plink2 --bfile "$BFILE" --pca 10 --out "${PREFIX}_pca"

echo -e "\n=== Step 2: Run Association ==="
plink2 --bfile "$BFILE" \
    --pheno "$PHENO" \
    --covar "${PREFIX}_pca.eigenvec" \
    --covar-name PC1-PC5 \
    --glm hide-covar \
    --out "$PREFIX"

echo -e "\n=== Step 3: Summarize Results ==="
RESULT_FILE=$(ls ${PREFIX}.*.glm.* 2>/dev/null | head -1)

if [[ -f "$RESULT_FILE" ]]; then
    echo "Results file: $RESULT_FILE"

    N_TESTS=$(tail -n +2 "$RESULT_FILE" | grep -c "ADD")
    N_GW=$(awk '$13 < 5e-8 && $7 == "ADD"' "$RESULT_FILE" | wc -l)
    N_SUG=$(awk '$13 < 1e-5 && $7 == "ADD"' "$RESULT_FILE" | wc -l)

    echo "Total tests: $N_TESTS"
    echo "Genome-wide significant (P < 5e-8): $N_GW"
    echo "Suggestive (P < 1e-5): $N_SUG"

    echo -e "\n=== Top 10 Hits ==="
    head -1 "$RESULT_FILE"
    awk '$7 == "ADD"' "$RESULT_FILE" | sort -k13 -g | head -10
else
    echo "No results file found"
    exit 1
fi

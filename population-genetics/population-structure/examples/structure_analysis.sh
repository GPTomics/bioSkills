#!/bin/bash
# Population structure analysis with PCA and Admixture
# Usage: ./structure_analysis.sh <plink_prefix> <output_prefix> [max_K]

BFILE="${1}"
PREFIX="${2:-structure}"
MAX_K="${3:-6}"

if [[ -z "$BFILE" ]]; then
    echo "Usage: $0 <plink_prefix> [output_prefix] [max_K]"
    exit 1
fi

echo "=== Population Structure Analysis ==="
echo "Input: $BFILE"
echo "Output prefix: $PREFIX"
echo "Testing K=2 to K=$MAX_K"

echo -e "\n=== Step 1: LD Pruning ==="
plink2 --bfile "$BFILE" --indep-pairwise 50 10 0.1 --out "${PREFIX}_prune"
plink2 --bfile "$BFILE" --extract "${PREFIX}_prune.prune.in" --make-bed --out "${PREFIX}_pruned"

echo "Pruned SNPs: $(wc -l < ${PREFIX}_pruned.bim)"

echo -e "\n=== Step 2: PCA ==="
plink2 --bfile "${PREFIX}_pruned" --pca 20 --out "${PREFIX}_pca"

echo -e "\n=== Step 3: Admixture ==="
cd "$(dirname ${PREFIX}_pruned.bed)" || exit

for K in $(seq 2 $MAX_K); do
    echo "Running K=$K..."
    admixture --cv -j4 "$(basename ${PREFIX}_pruned.bed)" $K 2>&1 | tee "${PREFIX}_log${K}.out"
done

echo -e "\n=== Step 4: CV Error Summary ==="
echo "K CV_Error"
for K in $(seq 2 $MAX_K); do
    CV=$(grep "CV" "${PREFIX}_log${K}.out" | awk '{print $4}')
    echo "$K $CV"
done

echo -e "\n=== Output Files ==="
echo "PCA eigenvectors: ${PREFIX}_pca.eigenvec"
echo "PCA eigenvalues: ${PREFIX}_pca.eigenval"
echo "Admixture Q files: ${PREFIX}_pruned.*.Q"

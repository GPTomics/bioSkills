#!/bin/bash
# LD analysis and pruning pipeline
# Usage: ./ld_analysis.sh <plink_prefix> <output_prefix>

BFILE="${1}"
PREFIX="${2:-ld_analysis}"

if [[ -z "$BFILE" ]]; then
    echo "Usage: $0 <plink_prefix> [output_prefix]"
    exit 1
fi

echo "=== LD Analysis Pipeline ==="
echo "Input: $BFILE"
echo "Output prefix: $PREFIX"

echo -e "\n=== Step 1: Calculate LD for QC ==="
plink2 --bfile "$BFILE" \
    --r2 \
    --ld-window-kb 500 \
    --ld-window-r2 0.2 \
    --out "${PREFIX}_ld"

N_PAIRS=$(wc -l < "${PREFIX}_ld.ld" 2>/dev/null || echo "0")
echo "LD pairs (r² > 0.2): $N_PAIRS"

echo -e "\n=== Step 2: Generate Pruned SNP Set ==="
plink2 --bfile "$BFILE" \
    --indep-pairwise 50 10 0.1 \
    --out "${PREFIX}_prune"

N_KEEP=$(wc -l < "${PREFIX}_prune.prune.in")
N_REMOVE=$(wc -l < "${PREFIX}_prune.prune.out")
echo "SNPs to keep: $N_KEEP"
echo "SNPs to remove: $N_REMOVE"

echo -e "\n=== Step 3: Create Pruned Dataset ==="
plink2 --bfile "$BFILE" \
    --extract "${PREFIX}_prune.prune.in" \
    --make-bed \
    --out "${PREFIX}_pruned"

echo -e "\n=== Step 4: Identify Haplotype Blocks ==="
plink --bfile "$BFILE" --blocks no-pheno-req --out "${PREFIX}_blocks" 2>/dev/null

if [[ -f "${PREFIX}_blocks.blocks.det" ]]; then
    N_BLOCKS=$(tail -n +2 "${PREFIX}_blocks.blocks.det" | wc -l)
    echo "Haplotype blocks identified: $N_BLOCKS"
fi

echo -e "\n=== Summary ==="
echo "Original SNPs: $(wc -l < ${BFILE}.bim)"
echo "Pruned SNPs: $(wc -l < ${PREFIX}_pruned.bim)"
echo "LD pairs file: ${PREFIX}_ld.ld"
echo "Pruned dataset: ${PREFIX}_pruned.{bed,bim,fam}"

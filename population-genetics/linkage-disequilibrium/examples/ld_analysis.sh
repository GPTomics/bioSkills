#!/usr/bin/env bash
# Reference: PLINK 1.9 (1.90b7+), PLINK 2.0 (alpha 6+) | Verify API if version differs
# LD pipeline that keeps the r2-vs-D' and pruning-vs-clumping boundaries explicit: exclude
# long-range-LD regions by COORDINATE, prune genotype-blind with composite r2, optionally clump
# GWAS summary stats p-value-aware, and define Gabriel blocks. Outputs go to a caller-supplied
# directory (default: a fresh temp dir) so nothing is written to the current dir.
# Usage: ./ld_analysis.sh <plink_prefix> [gwas_sumstats.txt] [output_dir]
set -euo pipefail

BFILE="${1:?Usage: $0 <plink_prefix> [gwas_sumstats.txt] [output_dir]}"
SUMSTATS="${2:-}"
OUTDIR="${3:-$(mktemp -d)}"
mkdir -p "$OUTDIR"
echo "Outputs -> $OUTDIR"

# Long-range-LD regions to special-case (hg19 coordinates; Price 2008). Their internal r2 is high
# and REAL, so a window prune cannot remove them - they must be excluded by position or they
# hijack the top PCs. Format: CHR START END LABEL.
RANGE="$OUTDIR/longrange_ld.txt"
cat > "$RANGE" <<'EOF'
6 25000000 34000000 MHC
8 8000000 12000000 inv8p23
17 40000000 45000000 inv17q21
EOF

# Step 1: drop long-range-LD regions by coordinate BEFORE pruning.
plink2 --bfile "$BFILE" --exclude range "$RANGE" --make-bed --out "$OUTDIR/noLR"

# Step 2: genotype-blind LD prune. 50-variant window, step 5, r2 0.1 (near-independence for PCA).
# Step must be 1 if the window were given in kb (plink2 rule); a variant-count window allows step 5.
plink2 --bfile "$OUTDIR/noLR" --indep-pairwise 50 5 0.1 --out "$OUTDIR/prune"
plink2 --bfile "$OUTDIR/noLR" --extract "$OUTDIR/prune.prune.in" --make-bed --out "$OUTDIR/pruned"
echo "Pruned-in variants: $(wc -l < "$OUTDIR/prune.prune.in")"

# Step 3: composite (unphased) r2 - robust default needing no phasing or HWE assumption.
plink2 --bfile "$BFILE" --r2-unphased --ld-window-kb 500 --ld-window-r2 0.2 --out "$OUTDIR/ld_unphased"

# Step 4 (optional): clump GWAS summary stats. p-value-aware, NOT a substitute for fine-mapping.
# Defaults (p1=1e-4, r2=0.5) are neither genome-wide nor strict, so they are overridden here.
if [[ -n "$SUMSTATS" ]]; then
    plink --bfile "$BFILE" --clump "$SUMSTATS" \
        --clump-p1 5e-8 --clump-p2 1e-5 --clump-r2 0.1 --clump-kb 250 \
        --out "$OUTDIR/clumped"
fi

# Step 5: Gabriel confidence-interval haplotype blocks (D'-based; for recombination/block maps).
plink --bfile "$BFILE" --blocks no-pheno-req --out "$OUTDIR/blocks"
if [[ -f "$OUTDIR/blocks.blocks.det" ]]; then
    echo "Haplotype blocks: $(($(wc -l < "$OUTDIR/blocks.blocks.det") - 1))"
fi

echo "Done. Pruned fileset: $OUTDIR/pruned.{bed,bim,fam}"

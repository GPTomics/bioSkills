#!/usr/bin/env bash
# Reference: PLINK 2.0 (alpha 6+), ADMIXTURE 1.3+ | Verify API if version differs
# Decision-grade population-structure pipeline: LD-prune + exclude long-range-LD/inversion regions,
# remove second-degree relatives BEFORE PCA, compute PCs, then run ADMIXTURE over a K SPAN with
# cross-validation (CV error is a GUIDE, never "the true K"). All outputs go to a caller-supplied
# directory (default: a fresh temp dir) so nothing is written to the current dir.
# Usage: ./structure_analysis.sh <plink_prefix> [max_K] [lrld_range_file] [output_dir]
set -euo pipefail

BFILE="${1:?Usage: $0 <plink_prefix> [max_K] [lrld_range_file] [output_dir]}"
MAX_K="${2:-6}"
LRLD="${3:-}"
OUTDIR="${4:-$(mktemp -d)}"
mkdir -p "$OUTDIR"
echo "Outputs -> $OUTDIR"

# LD-prune so a single dense block cannot dominate a PC (50-SNP window, 5-SNP step, r2 0.2).
plink2 --bfile "$BFILE" --indep-pairwise 50 5 0.2 --out "$OUTDIR/prune"

# Build the PCA input: extract pruned SNPs, drop very-low-MAF variants (destabilize PCA), and
# exclude long-range-LD/inversion regions (MHC, 8p23, 17q21.31, LCT) when a range file is supplied.
EXCLUDE=()
if [[ -n "$LRLD" ]]; then
    EXCLUDE=(--exclude range "$LRLD")
fi
plink2 --bfile "$BFILE" --extract "$OUTDIR/prune.prune.in" "${EXCLUDE[@]}" --maf 0.01 \
    --make-bed --out "$OUTDIR/for_pca"

# Remove relatives BEFORE computing axes: a cluster of relatives forms its own spurious PC.
# 0.0884 is the KING second-degree kinship cutoff (MZ 0.354 / 1st 0.177 / 2nd 0.0884 / 3rd 0.0442).
plink2 --bfile "$OUTDIR/for_pca" --king-cutoff 0.0884 --out "$OUTDIR/unrel"
plink2 --bfile "$OUTDIR/for_pca" --keep "$OUTDIR/unrel.king.cutoff.in.id" \
    --make-bed --out "$OUTDIR/pca_input"

# PCA on the pruned, inversion-stripped, unrelated set. approx is near-required above ~50k samples.
plink2 --bfile "$OUTDIR/pca_input" --pca 20 approx --out "$OUTDIR/pca"
echo "PCs -> $OUTDIR/pca.eigenvec  eigenvalues -> $OUTDIR/pca.eigenval"

# ADMIXTURE over a K SPAN with 5-fold cross-validation (the default). admixture writes .Q/.P next to
# its input, so run it inside OUTDIR on a basename to keep every artifact out of the current dir.
ADM_BED="$OUTDIR/pca_input.bed"
(
    cd "$OUTDIR"
    for K in $(seq 2 "$MAX_K"); do
        echo "ADMIXTURE K=$K ..."
        admixture --cv -j4 "$(basename "$ADM_BED")" "$K" 2>&1 | tee "log_K${K}.out"
    done
)

# CV error is a GUIDE: report the whole curve, not just the argmin. The true number of populations
# depends on sampling design and is often unidentifiable; present a span and check Q stability.
echo "=== Cross-validation error (guide, not the true K) ==="
grep -h "CV error" "$OUTDIR"/log_K*.out

echo "=== Output files ==="
echo "Eigenvectors: $OUTDIR/pca.eigenvec"
echo "Eigenvalues:  $OUTDIR/pca.eigenval"
echo "ADMIXTURE Q:  $OUTDIR/pca_input.*.Q"
echo "ADMIXTURE P:  $OUTDIR/pca_input.*.P"

#!/bin/bash
# Reference: ASTER 1.15+, IQ-TREE 2.2+ | Verify API if version differs
# Full ASTRAL species tree pipeline: gene tree inference, species tree estimation,
# and concordance factor computation.

set -euo pipefail

LOCI_DIR="loci"
OUTDIR="species_tree_results"
THREADS_PER_LOCUS=1
SEED=12345

mkdir -p "$OUTDIR"

# --- Step 1: Infer per-locus gene trees with IQ-TREE2 ---
# Each locus gets its own ML tree with model selection and ultrafast bootstrap.
# Single thread per locus allows running many loci in parallel via xargs/GNU parallel.
for f in "$LOCI_DIR"/*.fasta; do
    prefix="${OUTDIR}/$(basename "${f%.fasta}")"
    iqtree2 -s "$f" -m MFP -B 1000 -bnni -T "$THREADS_PER_LOCUS" \
        --prefix "$prefix" --seed "$SEED" --quiet
done

# --- Step 2: Collect all gene trees ---
cat "$OUTDIR"/*.treefile > "$OUTDIR/gene_trees.tre"
NGENES=$(wc -l < "$OUTDIR/gene_trees.tre")
echo "Collected $NGENES gene trees"

# --- Step 3: ASTRAL species tree estimation ---
# Local posterior probabilities annotated on branches (not bootstrap).
# Branch lengths are in coalescent units.
astral -i "$OUTDIR/gene_trees.tre" -o "$OUTDIR/species_tree.tre"
echo "Species tree: $OUTDIR/species_tree.tre"

# --- Step 4: wASTRAL (recommended when gene trees are noisy) ---
# Weights quartets by gene tree branch support, improving accuracy with
# short loci or low per-locus signal.
wastral -i "$OUTDIR/gene_trees.tre" -o "$OUTDIR/species_tree_wastral.tre"
echo "wASTRAL species tree: $OUTDIR/species_tree_wastral.tre"

# --- Step 5: Concordance factors ---
# Requires a concatenated alignment for site concordance factors.
# gCF: proportion of gene trees supporting each branch
# sCF: proportion of decisive sites supporting each branch (likelihood-based)
CONCAT="concat.fasta"
if [ -f "$CONCAT" ]; then
    iqtree2 -t "$OUTDIR/species_tree.tre" --gcf "$OUTDIR/gene_trees.tre" \
        -s "$CONCAT" --scfl 100 --prefix "$OUTDIR/concord"
    echo "Concordance factors: $OUTDIR/concord.cf.stat"
    echo "Annotated tree: $OUTDIR/concord.cf.tree"
else
    echo "Skipping concordance factors: $CONCAT not found"
    echo "Create concatenated alignment to compute gCF/sCF"
fi

echo "Pipeline complete."

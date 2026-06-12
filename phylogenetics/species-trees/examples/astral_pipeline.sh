#!/bin/bash
# Reference: ASTER 1.15+, IQ-TREE 2.2+ | Verify API if version differs
# Coalescent species-tree pipeline: per-locus gene trees -> contract weak branches ->
# wASTRAL/ASTRAL species tree -> gene and site concordance factors.
# Deliverable: a species tree consistent under the MSC PLUS the gCF/sCF that show where
# the data actually agree -- not a single bootstrap number that hides the discordance.
#
# NOT SPOT-RUNNABLE without ASTER (astral/wastral) + IQ-TREE2 installed and real loci.
# All outputs are written under a namespaced OUTDIR the caller can delete; nothing in CWD.

set -euo pipefail

LOCI_DIR="loci"
OUTDIR="species_tree_results"
CONTRACT_SUPPORT=10   # collapse gene-tree branches below 10% bootstrap to polytomies;
                      # widely-used default that cuts gene-tree-error bias (wASTRAL is the
                      # continuous-weighting alternative that removes this threshold choice)
THREADS=8
SEED=12345

mkdir -p "$OUTDIR"

# --- Step 1: per-locus gene trees with support (one ML tree per locus) ---
# -S runs IQ-TREE2 once per alignment in the directory; single concatenated treefile out.
iqtree2 -S "$LOCI_DIR" -m MFP -B 1000 -T AUTO --prefix "$OUTDIR/loci" --seed "$SEED" --quiet

# --- Step 2: contract weak gene-tree branches to polytomies ---
# ASTRAL-III handles polytomies correctly (they add no spurious quartet similarity).
# nw_ed (Newick Utilities): 'i & b<=N' selects internal nodes with support <= N to collapse.
if command -v nw_ed >/dev/null 2>&1; then
    nw_ed "$OUTDIR/loci.treefile" "i & b<=$CONTRACT_SUPPORT" o > "$OUTDIR/gene_trees.nwk"
else
    echo "nw_ed not found; using uncontracted gene trees (prefer wASTRAL instead)"
    cp "$OUTDIR/loci.treefile" "$OUTDIR/gene_trees.nwk"
fi
NGENES=$(wc -l < "$OUTDIR/gene_trees.nwk")
echo "Collected $NGENES gene trees"

# --- Step 3: primary estimate = wASTRAL (weights quartets by gene-tree support+length) ---
wastral -i "$OUTDIR/gene_trees.nwk" -o "$OUTDIR/species_wastral.tre" 2> "$OUTDIR/wastral.log"
echo "wASTRAL species tree: $OUTDIR/species_wastral.tre"

# --- Step 4: classic ASTRAL-III for localPP + full quartet annotation ---
# In ASTER, -t is THREADS and -u is annotation (the Java ASTRAL convention is the opposite).
# -u 2 annotates localPP and q1/q2/q3 for all three resolutions of each branch.
astral -t "$THREADS" -u 2 -i "$OUTDIR/gene_trees.nwk" -o "$OUTDIR/species_astral.tre" \
    2> "$OUTDIR/astral.log"
echo "ASTRAL species tree (localPP + q1/q2/q3): $OUTDIR/species_astral.tre"

# --- Step 5: gene and site concordance factors ---
# gCF = % of decisive gene trees containing each branch; sCF = % of decisive sites supporting it.
# A high-bootstrap branch with gCF ~ 25 is screaming disagreement the bootstrap hides.
CONCAT="concat.fasta"
if [ -f "$CONCAT" ]; then
    iqtree2 -te "$OUTDIR/species_astral.tre" --gcf "$OUTDIR/gene_trees.nwk" \
        -s "$CONCAT" --scfl 100 --prefix "$OUTDIR/cf"
    echo "Concordance factors: $OUTDIR/cf.cf.stat (per-branch q1/q2/q3)"
    echo "Annotated tree: $OUTDIR/cf.cf.tree"
else
    echo "Skipping concordance factors: $CONCAT not found"
    echo "Provide a concatenated alignment to compute gCF/sCF"
fi

echo "Pipeline complete. ASTRAL trees are UNROOTED; root with an outgroup, and remember"
echo "branch lengths are in coalescent units (not time). Trust the coalescent topology on low-gCF branches."

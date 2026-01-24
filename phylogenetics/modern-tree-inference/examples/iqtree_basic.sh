#!/bin/bash
# Basic IQ-TREE2 analysis with automatic model selection and ultrafast bootstrap

# Input: FASTA alignment
ALIGNMENT="alignment.fasta"
PREFIX="iqtree_analysis"

# For test data, download example alignment:
# wget https://raw.githubusercontent.com/Cibiv/IQ-TREE/master/example.phy

# Basic run with ModelFinder Plus and ultrafast bootstrap
# -m MFP: Automatic model selection integrated with tree inference
# -B 1000: 1000 UFBoot replicates (minimum for publication; use 10000 for final)
# -T AUTO: Detect available threads automatically
iqtree2 -s "$ALIGNMENT" -m MFP -B 1000 -T AUTO --prefix "$PREFIX"

# Output files:
# ${PREFIX}.treefile     - Best ML tree (Newick format)
# ${PREFIX}.contree      - Consensus tree with bootstrap support
# ${PREFIX}.iqtree       - Full report including model parameters
# ${PREFIX}.log          - Run log

echo "Best tree: ${PREFIX}.treefile"
echo "Report: ${PREFIX}.iqtree"

# View selected model
grep "Best-fit model" "${PREFIX}.iqtree"

# View tree
cat "${PREFIX}.treefile"

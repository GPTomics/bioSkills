#!/bin/bash
# Reference: IQ-TREE 2.2+ | Verify API if version differs
# Basic IQ-TREE2 analysis with automatic model selection and ultrafast bootstrap

# Input: FASTA alignment
ALIGNMENT="alignment.fasta"
PREFIX="iqtree_analysis"

# For test data, download example alignment:
# wget https://raw.githubusercontent.com/Cibiv/IQ-TREE/master/example.phy

# Standard analysis: ModelFinder + UFBoot2 + SH-aLRT
# -m MFP: ModelFinder Plus (tests FreeRate models; -m TEST does not)
# -B 1000: 1000 UFBoot replicates (minimum for publication; use 10000 for final)
# -alrt 1000: SH-aLRT for complementary support measure
# -bnni: Reduces UFBoot overestimation via NNI optimization of bootstrap trees
# --seed: Fixed seed for reproducibility
iqtree2 -s "$ALIGNMENT" -m MFP -B 1000 -alrt 1000 -bnni -T AUTO --seed 12345 --prefix "$PREFIX"

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

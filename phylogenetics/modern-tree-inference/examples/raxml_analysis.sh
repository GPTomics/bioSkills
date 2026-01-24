#!/bin/bash
# RAxML-ng maximum likelihood tree inference

ALIGNMENT="alignment.fasta"
PREFIX="raxml_analysis"

# Check alignment format and determine data type
raxml-ng --check --msa "$ALIGNMENT" --model GTR+G --prefix check

# Full analysis: ML search + bootstrap
# --all: Combined ML search and bootstrapping
# --model GTR+G: General Time Reversible + Gamma rate variation
# --bs-trees 100: 100 bootstrap replicates (use 1000 for publication)
# --threads auto: Automatic thread detection
raxml-ng --all \
    --msa "$ALIGNMENT" \
    --model GTR+G \
    --bs-trees 100 \
    --threads auto \
    --seed 12345 \
    --prefix "$PREFIX"

# Output files:
# ${PREFIX}.raxml.bestTree    - Best ML tree
# ${PREFIX}.raxml.support     - Tree with bootstrap support values
# ${PREFIX}.raxml.bootstraps  - Individual bootstrap trees
# ${PREFIX}.raxml.log         - Analysis log

echo "Best tree: ${PREFIX}.raxml.bestTree"
echo "With support: ${PREFIX}.raxml.support"

# For thorough ML search, use multiple starting trees
# --tree pars{10}: 10 parsimony starting trees
# --tree rand{10}: 10 random starting trees
raxml-ng --search \
    --msa "$ALIGNMENT" \
    --model GTR+G \
    --tree pars{10},rand{10} \
    --threads auto \
    --prefix thorough_search

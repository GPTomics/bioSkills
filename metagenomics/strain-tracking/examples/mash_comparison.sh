#!/bin/bash
# Strain comparison using MASH

set -euo pipefail

INPUT_DIR=${1:-.}
OUTPUT_PREFIX=${2:-mash_results}

echo "=== MASH Strain Comparison ==="

# Create sketches for all genomes
echo "Creating sketches..."
mash sketch -s 10000 -o ${OUTPUT_PREFIX}_sketch ${INPUT_DIR}/*.fasta

# Pairwise distances
echo "Calculating pairwise distances..."
mash dist ${OUTPUT_PREFIX}_sketch.msh ${OUTPUT_PREFIX}_sketch.msh > ${OUTPUT_PREFIX}_distances.tsv

# Find closely related strains (distance < 0.001)
echo ""
echo "=== Closely Related Strains (MASH distance < 0.001) ==="
awk '$1 != $2 && $3 < 0.001' ${OUTPUT_PREFIX}_distances.tsv | head -20

# Summary statistics
echo ""
echo "=== Distance Summary ==="
awk '$1 != $2 {print $3}' ${OUTPUT_PREFIX}_distances.tsv | sort -n | \
    awk '{a[NR]=$1} END {print "Min:", a[1]; print "Max:", a[NR]; print "Median:", a[int(NR/2)]}'

echo ""
echo "Results: ${OUTPUT_PREFIX}_distances.tsv"

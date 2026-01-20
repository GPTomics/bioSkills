#!/bin/bash
# AMR detection with AMRFinderPlus

set -euo pipefail

INPUT=${1:-contigs.fasta}
OUTPUT=${2:-amr_results.tsv}

echo "=== AMRFinderPlus AMR Detection ==="
echo "Input: $INPUT"

# Update database (run periodically)
# amrfinder -u

# Run AMRFinderPlus with extended output
amrfinder \
    -n ${INPUT} \
    -o ${OUTPUT} \
    --plus \
    --threads 8

# Summary statistics
echo ""
echo "=== Results Summary ==="
total=$(tail -n +2 ${OUTPUT} | wc -l)
echo "Total AMR genes found: $total"

echo ""
echo "By drug class:"
tail -n +2 ${OUTPUT} | cut -f12 | sort | uniq -c | sort -rn | head -10

echo ""
echo "By element type:"
tail -n +2 ${OUTPUT} | cut -f9 | sort | uniq -c | sort -rn

echo ""
echo "Results written to: ${OUTPUT}"

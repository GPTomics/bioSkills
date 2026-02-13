#!/bin/bash
# Reference: BioPython 1.83+, Infernal 1.1+, pandas 2.2+ | Verify API if version differs
# Search query sequences against the Rfam database using Infernal cmscan.

QUERY=$1
RFAM_CM=${2:-"Rfam.cm"}
RFAM_CLANIN=${3:-"Rfam.clanin"}
OUTPUT_PREFIX=${4:-"rfam_results"}
THREADS=${5:-8}

if [ -z "$QUERY" ]; then
    echo "Usage: $0 <query.fa> [Rfam.cm] [Rfam.clanin] [output_prefix] [threads]"
    echo ""
    echo "Searches query sequences against Rfam to classify ncRNA families."
    exit 1
fi

if [ ! -f "${RFAM_CM}.i1m" ]; then
    echo "Rfam CM not indexed. Running cmpress..."
    cmpress "$RFAM_CM"
fi

echo "=== Scanning $(grep -c '>' "$QUERY") sequences against Rfam ==="

cmscan \
    --cpu "$THREADS" \
    --tblout "${OUTPUT_PREFIX}.tbl" \
    --fmt 2 \
    --clanin "$RFAM_CLANIN" \
    --oclan \
    --cut_ga \
    --noali \
    "$RFAM_CM" \
    "$QUERY" > "${OUTPUT_PREFIX}.out"

echo "=== Results summary ==="
# Count significant hits (lines not starting with #)
NHITS=$(grep -cv '^#' "${OUTPUT_PREFIX}.tbl" 2>/dev/null || echo 0)
echo "Total hits: $NHITS"

echo ""
echo "Top families:"
grep -v '^#' "${OUTPUT_PREFIX}.tbl" | awk '{print $1}' | sort | uniq -c | sort -rn | head -10

echo ""
echo "Output files:"
echo "  ${OUTPUT_PREFIX}.tbl - Tabular results"
echo "  ${OUTPUT_PREFIX}.out - Full output with alignments"

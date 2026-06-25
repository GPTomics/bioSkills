#!/bin/bash
# Reference: Infernal 1.1.4+ | Verify API if version differs
# Classify ncRNA families by scanning sequences against Rfam with the documented pipeline.

QUERY=$1
RFAM_CM=${2:-"Rfam.cm"}
RFAM_CLANIN=${3:-"Rfam.clanin"}
OUTPUT_PREFIX=${4:-"rfam_results"}
THREADS=${5:-8}

if [ -z "$QUERY" ]; then
    echo "Usage: $0 <query.fa> [Rfam.cm] [Rfam.clanin] [output_prefix] [threads]"
    exit 1
fi

# Rfam.cm ships PRE-CALIBRATED: press it, never cmcalibrate it.
if [ ! -f "${RFAM_CM}.i1m" ]; then
    echo "Rfam CM not indexed; running cmpress..."
    cmpress "$RFAM_CM"
fi

# --cut_ga uses the curated per-family gathering thresholds (bit scores, DB-size-independent),
# the correct default over a flat -E. --rfam is the large-DB strict filter; --nohmmonly keeps GA
# valid for every model. -Z <2 x total Mb> makes E-values reproducible (compute from the query).
DBMB=$(awk '!/^>/{n+=length($0)} END{printf "%.4f", (n*2)/1000000}' "$QUERY")
echo "=== Scanning $(grep -c '>' "$QUERY") sequences against Rfam (-Z ${DBMB} Mb) ==="

cmscan \
    --cpu "$THREADS" \
    -Z "$DBMB" \
    --cut_ga \
    --rfam \
    --nohmmonly \
    --fmt 2 \
    --clanin "$RFAM_CLANIN" \
    --tblout "${OUTPUT_PREFIX}.tbl" \
    "$RFAM_CM" \
    "$QUERY" > "${OUTPUT_PREFIX}.out"

# Clan deoverlapping: drop hits marked '=' (dominated by a higher-scoring clanmate).
grep -v ' = ' "${OUTPUT_PREFIX}.tbl" > "${OUTPUT_PREFIX}.deoverlapped.tbl"

echo "=== Results ==="
NHITS=$(grep -cv '^#' "${OUTPUT_PREFIX}.deoverlapped.tbl" 2>/dev/null || echo 0)
echo "Hits after deoverlapping: $NHITS"
echo "Top families:"
grep -v '^#' "${OUTPUT_PREFIX}.deoverlapped.tbl" | awk '{print $2}' | sort | uniq -c | sort -rn | head -10

echo ""
echo "Output: ${OUTPUT_PREFIX}.deoverlapped.tbl (parse with parse_infernal.py)"

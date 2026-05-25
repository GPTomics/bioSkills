#!/bin/bash
# Reference: NCBI Datasets CLI 16.0+, aria2c 1.36+ | Verify API if version differs
# Bulk pull via --dehydrated + parallel transfer (the cloud / HPC pattern).

set -euo pipefail

TAXON="${1:-Salmonella enterica}"
OUT_ZIP="${2:-bulk_dehydrated.zip}"
DEST="${3:-bulk_dataset}"
THREADS="${4:-8}"

echo "=== Step 1: dehydrated discovery (small ZIP, just metadata) ==="
datasets download genome taxon "${TAXON}" \
    --reference \
    --annotated \
    --assembly-source RefSeq \
    --include genome,gff3,protein \
    --dehydrated \
    --filename "${OUT_ZIP}" \
    --no-progressbar

unzip -q "${OUT_ZIP}" -d "${DEST}/"

FETCH="${DEST}/ncbi_dataset/fetch.txt"
echo "  Files queued: $(wc -l < ${FETCH})"

echo
echo "=== Step 2: parallel transfer with aria2c ==="
# The fetch.txt format is: <url> [TAB] <path-relative-to-data-dir>
# Older aria2c doesn't accept this format; transform to input file:
awk -F'\t' '{print $1"\n  out="$2}' "${FETCH}" > "${DEST}/aria2_input.txt"

aria2c \
    --input-file="${DEST}/aria2_input.txt" \
    --dir="${DEST}/ncbi_dataset/data/" \
    --max-concurrent-downloads="${THREADS}" \
    --max-connection-per-server="${THREADS}" \
    --retry-wait=5 \
    --quiet=true

echo
echo "=== Step 3: verify with datasets rehydrate ==="
# datasets rehydrate validates checksums of all files
datasets rehydrate --directory "${DEST}/" --max-workers "${THREADS}"

echo
echo "=== Done ==="
ls "${DEST}/ncbi_dataset/data/" | head
echo
echo "Files: $(find ${DEST}/ncbi_dataset/data/ -type f | wc -l)"
echo "Total size: $(du -sh ${DEST}/ncbi_dataset/data/ | cut -f1)"

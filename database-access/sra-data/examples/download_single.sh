#!/bin/bash
# Reference: sra-tools 3.0+ | Verify API if version differs
# Single-run SRA download via the toolkit path: prefetch (with explicit --max-size) + vdb-validate + fasterq-dump + pigz.

set -euo pipefail

SRR="${1:-SRR12345678}"
OUT="${2:-./fastq}"
THREADS="${3:-8}"
MAX_SIZE="${4:-100G}"  # Default 20G silently skips larger -- always set explicitly

mkdir -p "${OUT}"

echo "=== prefetch ${SRR} (max-size ${MAX_SIZE}) ==="
prefetch "${SRR}" --max-size "${MAX_SIZE}" -p

echo
echo "=== vdb-validate ==="
vdb-validate "${SRR}" || { echo "VALIDATION FAILED"; exit 1; }

echo
echo "=== fasterq-dump (writes uncompressed; needs ~3x final size in scratch) ==="
# --split-files: emit _1.fastq and _2.fastq for paired
# DROP --skip-technical if this is 10x or other single-cell data (need barcodes/UMIs)
fasterq-dump "${SRR}" \
    -O "${OUT}" \
    -e "${THREADS}" \
    -p \
    --split-files \
    --skip-technical

echo
echo "=== pigz compression (fasterq-dump does NOT compress) ==="
pigz -p "${THREADS}" "${OUT}/${SRR}"_*.fastq

echo
echo "Files:"
ls -lh "${OUT}/${SRR}"_*.fastq.gz

echo
echo "Optional cleanup:"
echo "  rm -rf \"${SRR}\"   # remove cached .sra after successful FASTQ extraction"

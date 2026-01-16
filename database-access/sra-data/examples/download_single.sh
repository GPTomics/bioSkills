#!/bin/bash
# Download a single SRA run as FASTQ

SRR="${1:-SRR12345678}"
OUTDIR="${2:-./fastq}"
THREADS="${3:-4}"

echo "Downloading $SRR to $OUTDIR with $THREADS threads"

mkdir -p "$OUTDIR"

fasterq-dump "$SRR" \
    -O "$OUTDIR" \
    -e "$THREADS" \
    -p \
    --skip-technical

echo "Done. Files:"
ls -lh "${OUTDIR}"/${SRR}*

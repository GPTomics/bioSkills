#!/bin/bash
# Prefetch large SRA files before converting (more reliable)

SRR="${1:-SRR12345678}"
OUTDIR="${2:-./fastq}"
THREADS="${3:-8}"

echo "Prefetching $SRR..."
prefetch "$SRR" -p

echo "Validating..."
if vdb-validate "$SRR"; then
    echo "Validation passed"
else
    echo "Validation failed!"
    exit 1
fi

echo "Converting to FASTQ..."
mkdir -p "$OUTDIR"
fasterq-dump "$SRR" -O "$OUTDIR" -e "$THREADS" -p

echo "Cleaning up .sra file..."
rm -rf ~/ncbi/sra/${SRR}*

echo "Done. Files in $OUTDIR:"
ls -lh "${OUTDIR}"/${SRR}*

#!/bin/bash
# Long-read assembly with Flye
set -euo pipefail

READS=$1
OUTDIR=$2
SIZE=${3:-5m}
THREADS=${4:-16}

echo "=== Long-Read Assembly with Flye ==="
echo "Reads: $READS"
echo "Genome size: $SIZE"

# Determine read type from filename
if [[ "$READS" == *"hifi"* ]] || [[ "$READS" == *"ccs"* ]]; then
    READ_TYPE="--pacbio-hifi"
elif [[ "$READS" == *"pacbio"* ]] || [[ "$READS" == *"pb"* ]]; then
    READ_TYPE="--pacbio-raw"
else
    READ_TYPE="--nano-raw"
fi

echo "Read type: $READ_TYPE"

# Run Flye
flye \
    $READ_TYPE $READS \
    --out-dir $OUTDIR \
    --genome-size $SIZE \
    --threads $THREADS

# Output stats
echo ""
echo "=== Assembly Statistics ==="
cat ${OUTDIR}/assembly_info.txt

echo ""
echo "Assembly complete: ${OUTDIR}/assembly.fasta"

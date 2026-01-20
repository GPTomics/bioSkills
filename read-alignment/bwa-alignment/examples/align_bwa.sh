#!/bin/bash
# DNA alignment with BWA-MEM2

set -euo pipefail

R1=${1:-reads_1.fastq.gz}
R2=${2:-reads_2.fastq.gz}
REFERENCE=${3:-reference.fa}
OUTPUT=${4:-aligned.bam}
SAMPLE=${5:-sample}
THREADS=${6:-8}

echo "=== BWA-MEM2 Alignment ==="

# Index reference if needed
if [[ ! -f "${REFERENCE}.bwt.2bit.64" ]]; then
    echo "Indexing reference..."
    bwa-mem2 index $REFERENCE
fi

# Align with read group
bwa-mem2 mem -t $THREADS \
    -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" \
    $REFERENCE $R1 $R2 | \
    samtools sort -@ $THREADS -o $OUTPUT -

# Index BAM
samtools index $OUTPUT

echo "Output: $OUTPUT"
echo "$(samtools view -c $OUTPUT) reads aligned"

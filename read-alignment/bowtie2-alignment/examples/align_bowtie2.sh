#!/bin/bash
# Alignment with Bowtie2 (ATAC-seq, ChIP-seq)

set -euo pipefail

R1=${1:-reads_1.fastq.gz}
R2=${2:-reads_2.fastq.gz}
INDEX=${3:-bt2_index}
OUTPUT=${4:-aligned.bam}
THREADS=${5:-8}

echo "=== Bowtie2 Alignment ==="

# Align with parameters suitable for ATAC/ChIP-seq
bowtie2 -p $THREADS \
    -x $INDEX \
    -1 $R1 -2 $R2 \
    --very-sensitive \
    --no-mixed \
    --no-discordant \
    -X 2000 \
    2> ${OUTPUT%.bam}_stats.txt | \
    samtools view -@ $THREADS -bS - | \
    samtools sort -@ $THREADS -o $OUTPUT -

samtools index $OUTPUT

echo "Output: $OUTPUT"
cat ${OUTPUT%.bam}_stats.txt

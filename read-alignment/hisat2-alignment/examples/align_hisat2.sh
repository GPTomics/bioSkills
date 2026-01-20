#!/bin/bash
# RNA-seq alignment with HISAT2

set -euo pipefail

R1=${1:-reads_1.fastq.gz}
R2=${2:-reads_2.fastq.gz}
INDEX=${3:-hisat2_index}
OUTPUT=${4:-aligned.bam}
THREADS=${5:-8}

echo "=== HISAT2 RNA-seq Alignment ==="

hisat2 -p $THREADS \
    -x $INDEX \
    -1 $R1 -2 $R2 \
    --dta \
    --new-summary \
    --summary-file ${OUTPUT%.bam}_summary.txt \
    2> ${OUTPUT%.bam}_log.txt | \
    samtools sort -@ $THREADS -o $OUTPUT -

samtools index $OUTPUT

echo "Output: $OUTPUT"
cat ${OUTPUT%.bam}_summary.txt

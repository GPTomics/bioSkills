#!/bin/bash

REFERENCE="reference.fa"
READS="reads.fastq.gz"
OUTPUT="aligned"
THREADS=8
SAMPLE="sample"

minimap2 -d ${REFERENCE}.mmi $REFERENCE

minimap2 -ax map-ont \
    -t $THREADS \
    -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}" \
    --MD \
    ${REFERENCE}.mmi $READS | \
    samtools sort -@ 4 -o ${OUTPUT}.bam

samtools index ${OUTPUT}.bam

echo "Alignment statistics:"
samtools flagstat ${OUTPUT}.bam

echo ""
echo "Coverage estimate:"
samtools depth ${OUTPUT}.bam | awk '{sum+=$3} END {print "Mean depth: " sum/NR "x"}'

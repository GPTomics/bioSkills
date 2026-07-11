#!/bin/bash
# Reference: bwa 0.7.17+, samtools 1.19+ | Verify API if version differs
# Align, sort, and index reads

set -e

REF=$1
R1=$2
OUTPUT=$3
R2=$4
THREADS=${5:-8}

if [ -z "$REF" ] || [ -z "$R1" ] || [ -z "$OUTPUT" ]; then
    echo "Usage: sort_pipeline.sh <reference.fa> <R1.fq> <output.bam> [R2.fq] [threads]"
    echo "Example (paired): sort_pipeline.sh ref.fa reads_R1.fq aligned.bam reads_R2.fq 8"
    echo "Example (single): sort_pipeline.sh ref.fa reads.fq aligned.bam"
    exit 1
fi

echo "Aligning and sorting..."

if [ -z "$R2" ]; then
    # Single-end (R2 omitted)
    bwa mem -t "$THREADS" "$REF" "$R1" | samtools sort -@ 4 -o "$OUTPUT"
else
    # Paired-end
    bwa mem -t "$THREADS" "$REF" "$R1" "$R2" | samtools sort -@ 4 -o "$OUTPUT"
fi

echo "Indexing..."
samtools index "$OUTPUT"

echo "Done: $OUTPUT"
samtools flagstat "$OUTPUT"

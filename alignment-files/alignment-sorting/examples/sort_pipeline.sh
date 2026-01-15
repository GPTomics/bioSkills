#!/bin/bash
# Align, sort, and index reads

set -e

REF=$1
R1=$2
R2=$3
OUTPUT=$4
THREADS=${5:-8}

if [ -z "$REF" ] || [ -z "$R1" ] || [ -z "$OUTPUT" ]; then
    echo "Usage: sort_pipeline.sh <reference.fa> <R1.fq> [R2.fq] <output.bam> [threads]"
    echo "Example: sort_pipeline.sh ref.fa reads_R1.fq reads_R2.fq aligned.bam 8"
    exit 1
fi

echo "Aligning and sorting..."

if [ -z "$R2" ] || [ "$R2" == "$OUTPUT" ]; then
    # Single-end
    bwa mem -t "$THREADS" "$REF" "$R1" | samtools sort -@ 4 -o "$OUTPUT"
else
    # Paired-end
    bwa mem -t "$THREADS" "$REF" "$R1" "$R2" | samtools sort -@ 4 -o "$OUTPUT"
fi

echo "Indexing..."
samtools index "$OUTPUT"

echo "Done: $OUTPUT"
samtools flagstat "$OUTPUT"

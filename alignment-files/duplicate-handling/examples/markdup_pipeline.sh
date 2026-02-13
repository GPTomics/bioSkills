#!/bin/bash
# Reference: picard 3.1+, pysam 0.22+, samtools 1.19+ | Verify API if version differs
# Mark duplicates using samtools pipeline

set -e

INPUT=$1
OUTPUT=$2
THREADS=${3:-4}

if [ -z "$INPUT" ] || [ -z "$OUTPUT" ]; then
    echo "Usage: markdup_pipeline.sh <input.bam> <output.bam> [threads]"
    exit 1
fi

echo "Marking duplicates in $INPUT..."
echo "Using $THREADS threads"

samtools sort -n -@ "$THREADS" "$INPUT" | \
    samtools fixmate -m -@ "$THREADS" - - | \
    samtools sort -@ "$THREADS" - | \
    samtools markdup -@ "$THREADS" -s - "$OUTPUT" 2>&1 | tee markdup_stats.txt

echo "Indexing..."
samtools index "$OUTPUT"

echo "Done: $OUTPUT"
echo ""
echo "Duplicate statistics:"
samtools flagstat "$OUTPUT" | grep -E "(total|duplicates)"

#!/bin/bash
# Convert between SAM/BAM/CRAM formats

set -e

INPUT=$1
OUTPUT=$2
REFERENCE=$3

if [ -z "$INPUT" ] || [ -z "$OUTPUT" ]; then
    echo "Usage: convert_formats.sh <input> <output> [reference.fa]"
    echo "Examples:"
    echo "  convert_formats.sh input.sam output.bam"
    echo "  convert_formats.sh input.bam output.cram reference.fa"
    exit 1
fi

EXT="${OUTPUT##*.}"

case "$EXT" in
    sam)
        samtools view -h -o "$OUTPUT" "$INPUT"
        ;;
    bam)
        samtools view -b -o "$OUTPUT" "$INPUT"
        ;;
    cram)
        if [ -z "$REFERENCE" ]; then
            echo "Error: CRAM conversion requires reference.fa"
            exit 1
        fi
        samtools view -C -T "$REFERENCE" -o "$OUTPUT" "$INPUT"
        ;;
    *)
        echo "Unknown output format: $EXT"
        exit 1
        ;;
esac

echo "Converted $INPUT -> $OUTPUT"

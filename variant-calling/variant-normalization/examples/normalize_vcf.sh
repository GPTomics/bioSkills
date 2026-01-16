#!/bin/bash
# Normalize VCF file: left-align indels and split multiallelic sites

set -e

if [ $# -lt 3 ]; then
    echo "Usage: $0 <reference.fa> <input.vcf.gz> <output.vcf.gz>"
    exit 1
fi

REF="$1"
INPUT="$2"
OUTPUT="$3"

# Check inputs
if [ ! -f "$REF" ]; then
    echo "Error: Reference not found: $REF"
    exit 1
fi

if [ ! -f "$INPUT" ]; then
    echo "Error: Input VCF not found: $INPUT"
    exit 1
fi

# Count before
BEFORE=$(bcftools view -H "$INPUT" | wc -l | tr -d ' ')

# Normalize
echo "Normalizing VCF..."
bcftools norm -f "$REF" -m-any -d exact "$INPUT" -Oz -o "$OUTPUT"
bcftools index "$OUTPUT"

# Count after
AFTER=$(bcftools view -H "$OUTPUT" | wc -l | tr -d ' ')

# Report
echo ""
echo "=== Normalization Complete ==="
echo "Input:  $INPUT"
echo "Output: $OUTPUT"
echo ""
echo "Variants before: $BEFORE"
echo "Variants after:  $AFTER"
echo "Difference:      $((AFTER - BEFORE))"

#!/bin/bash
# Generate consensus sequence from VCF

set -e

if [ $# -lt 3 ]; then
    echo "Usage: $0 <reference.fa> <input.vcf.gz> <output.fa> [sample]"
    exit 1
fi

REF="$1"
VCF="$2"
OUTPUT="$3"
SAMPLE="${4:-}"

# Check inputs
if [ ! -f "$REF" ]; then
    echo "Error: Reference not found: $REF"
    exit 1
fi

if [ ! -f "$VCF" ]; then
    echo "Error: VCF not found: $VCF"
    exit 1
fi

# Check VCF is indexed
if [ ! -f "${VCF}.csi" ] && [ ! -f "${VCF}.tbi" ]; then
    echo "Indexing VCF..."
    bcftools index "$VCF"
fi

# Build command
CMD="bcftools consensus -f $REF"

if [ -n "$SAMPLE" ]; then
    CMD="$CMD -s $SAMPLE"
fi

CMD="$CMD $VCF"

# Generate consensus
echo "Generating consensus..."
$CMD > "$OUTPUT"

# Report
echo ""
echo "=== Consensus Generated ==="
echo "Reference: $REF"
echo "VCF: $VCF"
if [ -n "$SAMPLE" ]; then
    echo "Sample: $SAMPLE"
fi
echo "Output: $OUTPUT"
echo ""
echo "Variants applied: $(bcftools view -H "$VCF" | wc -l | tr -d ' ')"

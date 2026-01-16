#!/bin/bash
# Annotate VCF with rsIDs and population frequencies

set -e

if [ $# -lt 3 ]; then
    echo "Usage: $0 <input.vcf.gz> <dbsnp.vcf.gz> <output.vcf.gz>"
    echo ""
    echo "Optional: Set GNOMAD_VCF environment variable to add population frequencies"
    exit 1
fi

INPUT="$1"
DBSNP="$2"
OUTPUT="$3"

# Check inputs
for vcf in "$INPUT" "$DBSNP"; do
    if [ ! -f "$vcf" ]; then
        echo "Error: File not found: $vcf"
        exit 1
    fi
done

echo "Adding rsIDs from dbSNP..."
if [ -n "$GNOMAD_VCF" ] && [ -f "$GNOMAD_VCF" ]; then
    echo "Adding population frequencies from gnomAD..."
    bcftools annotate -a "$DBSNP" -c ID "$INPUT" | \
        bcftools annotate -a "$GNOMAD_VCF" -c INFO/AF -Oz -o "$OUTPUT"
else
    bcftools annotate -a "$DBSNP" -c ID "$INPUT" -Oz -o "$OUTPUT"
fi

bcftools index "$OUTPUT"

# Report
TOTAL=$(bcftools view -H "$OUTPUT" | wc -l | tr -d ' ')
WITH_ID=$(bcftools view -H "$OUTPUT" | awk -F'\t' '$3!="."' | wc -l | tr -d ' ')

echo ""
echo "=== Annotation Complete ==="
echo "Output: $OUTPUT"
echo "Total variants: $TOTAL"
echo "With rsID: $WITH_ID ($(echo "scale=1; $WITH_ID * 100 / $TOTAL" | bc)%)"

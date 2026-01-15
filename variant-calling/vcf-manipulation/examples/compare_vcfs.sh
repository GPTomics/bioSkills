#!/bin/bash
# Compare two VCF files and report overlap statistics

set -e

if [ $# -lt 2 ]; then
    echo "Usage: $0 <vcf1.vcf.gz> <vcf2.vcf.gz> [output_dir]"
    exit 1
fi

VCF1="$1"
VCF2="$2"
OUTDIR="${3:-vcf_comparison}"

# Check inputs
for vcf in "$VCF1" "$VCF2"; do
    if [ ! -f "$vcf" ]; then
        echo "Error: File not found: $vcf"
        exit 1
    fi
done

# Create output directory
mkdir -p "$OUTDIR"

# Run intersection
echo "Comparing VCF files..."
bcftools isec -p "$OUTDIR" -Oz "$VCF1" "$VCF2"

# Count variants in each category
VCF1_ONLY=$(bcftools view -H "${OUTDIR}/0000.vcf.gz" | wc -l | tr -d ' ')
VCF2_ONLY=$(bcftools view -H "${OUTDIR}/0001.vcf.gz" | wc -l | tr -d ' ')
SHARED=$(bcftools view -H "${OUTDIR}/0002.vcf.gz" | wc -l | tr -d ' ')

TOTAL1=$((VCF1_ONLY + SHARED))
TOTAL2=$((VCF2_ONLY + SHARED))

# Report
echo ""
echo "=== Comparison Results ==="
echo "VCF1: $VCF1"
echo "VCF2: $VCF2"
echo ""
echo "VCF1 only:  $VCF1_ONLY"
echo "VCF2 only:  $VCF2_ONLY"
echo "Shared:     $SHARED"
echo ""
echo "VCF1 total: $TOTAL1"
echo "VCF2 total: $TOTAL2"
echo ""
if [ "$TOTAL1" -gt 0 ]; then
    echo "Overlap with VCF1: $(echo "scale=1; $SHARED * 100 / $TOTAL1" | bc)%"
fi
if [ "$TOTAL2" -gt 0 ]; then
    echo "Overlap with VCF2: $(echo "scale=1; $SHARED * 100 / $TOTAL2" | bc)%"
fi
echo ""
echo "Output files in: $OUTDIR"

#!/bin/bash
# Reference: bcftools 1.12+ | Verify API if version differs
# Normalize a VCF (atomize MNPs -> split multiallelic -> left-align + parsimony)
# and quantify how many extra records atomization alone contributes, since
# vt decompose_blocksub splits MNPs by default while bcftools norm does not: mixing
# tools across cohorts manufactures spurious cohort-private variants at every MNP.

set -euo pipefail

if [ $# -lt 3 ]; then
    echo "Usage: $0 <reference.fa> <input.vcf.gz> <output.vcf.gz>"
    exit 1
fi

REF="$1"
INPUT="$2"
OUTPUT="$3"

if [ ! -f "$REF" ]; then
    echo "Error: Reference not found: $REF"
    exit 1
fi

if [ ! -f "$INPUT" ]; then
    echo "Error: Input VCF not found: $INPUT"
    exit 1
fi

count_records() { bcftools view -H "$1" | wc -l | tr -d ' '; }

BEFORE=$(count_records "$INPUT")

# Split + left-align WITHOUT atomization: MNPs stay as single records.
NO_ATOMIZE=$(bcftools norm -m- -f "$REF" "$INPUT" -Ou 2>/dev/null | bcftools view -H | wc -l | tr -d ' ')

# Full canonical pipeline: atomize MNPs, split multiallelic, left-align + parsimony.
bcftools norm --atomize "$INPUT" -Ou 2>/dev/null \
    | bcftools norm -m- -Ou \
    | bcftools norm -f "$REF" -Oz -o "$OUTPUT"
bcftools index -f "$OUTPUT"

AFTER=$(count_records "$OUTPUT")

echo ""
echo "=== Normalization Complete ==="
echo "Input:  $INPUT"
echo "Output: $OUTPUT"
echo ""
echo "Records before:                    $BEFORE"
echo "Records after split+left-align:    $NO_ATOMIZE"
echo "Records after full (with atomize): $AFTER"
echo "Extra records from atomization:    $((AFTER - NO_ATOMIZE))"
echo ""
echo "The difference counts the MNP and complex records that bcftools-only norm leaves"
echo "intact but vt-style decomposition splits into atoms. Standardize one tool +"
echo "flags across every cohort compared, or these become false private variants."

#!/bin/bash
# Reference: bcftools 1.19+, samtools 1.19+, bedtools 2.31+ | Verify API if version differs
# Generate a consensus FASTA by applying VCF variants onto a reference, with
# normalization and optional no-coverage masking built from callable depth.

set -euo pipefail

if [ $# -lt 3 ]; then
    echo "Usage: $0 <reference.fa> <input.vcf.gz> <output.fa> [sample] [bam] [min_depth]"
    echo "  bam + min_depth (optional): mask positions below min_depth to N so"
    echo "  no-coverage sites are not silently emitted as reference."
    exit 1
fi

REF="$1"; VCF="$2"; OUTPUT="$3"; SAMPLE="${4:-}"; BAM="${5:-}"; MIN_DEPTH="${6:-10}"

[ -f "$REF" ] || { echo "Error: reference not found: $REF"; exit 1; }
[ -f "$VCF" ] || { echo "Error: VCF not found: $VCF"; exit 1; }

# bcftools consensus requires a bgzipped, indexed VCF.
if [ ! -f "${VCF}.csi" ] && [ ! -f "${VCF}.tbi" ]; then
    echo "Indexing VCF..."
    bcftools index "$VCF"
fi

# Normalize first: left-align indels and split multiallelics so records match the
# reference context. Un-normalized indels produce wrong sequence with only a stderr warning.
NORM=$(mktemp -u).vcf.gz
bcftools norm -f "$REF" "$VCF" -Oz -o "$NORM"
bcftools index "$NORM"

# Build the consensus command.
CMD=(bcftools consensus -f "$REF")
[ -n "$SAMPLE" ] && CMD+=(-s "$SAMPLE")

# Optional mask: samtools depth -a is mandatory -- without -a, zero-coverage
# positions are OMITTED from the output and would escape masking, staying as
# reference. min_depth default 10x is a common floor for a confident base call.
if [ -n "$BAM" ]; then
    MASK=$(mktemp).bed
    samtools depth -a "$BAM" | awk -v d="$MIN_DEPTH" '$3 < d {print $1"\t"$2-1"\t"$2}' | bedtools merge > "$MASK"
    CMD+=(-m "$MASK")
    echo "Masking positions below ${MIN_DEPTH}x from $BAM"
fi

echo "Generating consensus..."
"${CMD[@]}" "$NORM" > "$OUTPUT"

echo ""
echo "=== Consensus generated ==="
echo "Reference: $REF | VCF: $VCF${SAMPLE:+ | Sample: $SAMPLE} | Output: $OUTPUT"
echo "Variants available to apply: $(bcftools view -H "$NORM" | wc -l | tr -d ' ')"

rm -f "$NORM" "${NORM}.csi" "${MASK:-}"

#!/bin/bash
# Call variants from BAM file using bcftools

set -e

if [ $# -lt 3 ]; then
    echo "Usage: $0 <reference.fa> <input.bam> <output.vcf.gz>"
    exit 1
fi

REF="$1"
BAM="$2"
OUTPUT="$3"

# Check inputs exist
if [ ! -f "$REF" ]; then
    echo "Error: Reference not found: $REF"
    exit 1
fi

if [ ! -f "$BAM" ]; then
    echo "Error: BAM not found: $BAM"
    exit 1
fi

# Check reference is indexed
if [ ! -f "${REF}.fai" ]; then
    echo "Indexing reference..."
    samtools faidx "$REF"
fi

# Call variants
echo "Calling variants..."
bcftools mpileup -Ou -f "$REF" \
    -q 20 -Q 20 \
    -a FORMAT/DP,FORMAT/AD \
    "$BAM" | \
bcftools call -mv -Oz -o "$OUTPUT"

# Index output
echo "Indexing output..."
bcftools index "$OUTPUT"

# Print summary
echo ""
echo "Variant calling complete: $OUTPUT"
echo ""
bcftools stats "$OUTPUT" | grep -E "^SN|number of"

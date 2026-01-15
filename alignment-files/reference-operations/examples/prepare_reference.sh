#!/bin/bash
# Prepare reference genome for analysis

set -e

REF=$1

if [ -z "$REF" ]; then
    echo "Usage: prepare_reference.sh <reference.fa>"
    exit 1
fi

if [ ! -f "$REF" ]; then
    echo "Error: Reference file not found: $REF"
    exit 1
fi

echo "Preparing reference: $REF"

echo "1. Creating FASTA index..."
samtools faidx "$REF"

echo "2. Creating sequence dictionary..."
samtools dict "$REF" -o "${REF%.fa}.dict"

echo "3. Generating chromosome sizes..."
cut -f1,2 "${REF}.fai" > "${REF%.fa}.chrom.sizes"

echo ""
echo "Files created:"
ls -la "${REF}"* "${REF%.fa}".*

echo ""
echo "Reference summary:"
echo "  Chromosomes: $(wc -l < "${REF}.fai")"
echo "  Total length: $(awk '{sum += $2} END {print sum}' "${REF}.fai") bp"

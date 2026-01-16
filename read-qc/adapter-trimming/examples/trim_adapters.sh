#!/bin/bash
# Trim Illumina adapters from paired-end reads
# Usage: ./trim_adapters.sh <R1.fastq.gz> <R2.fastq.gz> <output_prefix>

R1="${1}"
R2="${2}"
PREFIX="${3:-trimmed}"
THREADS="${4:-4}"

if [[ -z "$R1" ]] || [[ -z "$R2" ]]; then
    echo "Usage: $0 <R1.fastq.gz> <R2.fastq.gz> [output_prefix] [threads]"
    exit 1
fi

echo "=== Adapter Trimming ==="
echo "Read 1: $R1"
echo "Read 2: $R2"
echo "Output prefix: $PREFIX"

ADAPT_R1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPT_R2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

cutadapt \
    -a "$ADAPT_R1" \
    -A "$ADAPT_R2" \
    -m 20 \
    -j "$THREADS" \
    --report=minimal \
    -o "${PREFIX}_R1.fastq.gz" \
    -p "${PREFIX}_R2.fastq.gz" \
    "$R1" "$R2"

if [[ $? -eq 0 ]]; then
    echo -e "\n=== Output ==="
    echo "Trimmed R1: ${PREFIX}_R1.fastq.gz"
    echo "Trimmed R2: ${PREFIX}_R2.fastq.gz"

    echo -e "\n=== Read Counts ==="
    echo "Input R1: $(zcat "$R1" | wc -l | awk '{print $1/4}')"
    echo "Output R1: $(zcat "${PREFIX}_R1.fastq.gz" | wc -l | awk '{print $1/4}')"
else
    echo "Trimming failed"
    exit 1
fi

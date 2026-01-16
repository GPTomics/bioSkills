#!/bin/bash
# Quality filter paired-end reads with Trimmomatic
# Usage: ./quality_filter.sh <R1.fastq.gz> <R2.fastq.gz> <output_prefix>

R1="${1}"
R2="${2}"
PREFIX="${3:-filtered}"
THREADS="${4:-4}"

if [[ -z "$R1" ]] || [[ -z "$R2" ]]; then
    echo "Usage: $0 <R1.fastq.gz> <R2.fastq.gz> [output_prefix] [threads]"
    exit 1
fi

echo "=== Quality Filtering ==="
echo "Input R1: $R1"
echo "Input R2: $R2"
echo "Output prefix: $PREFIX"

trimmomatic PE -threads "$THREADS" -phred33 \
    "$R1" "$R2" \
    "${PREFIX}_R1.fastq.gz" "${PREFIX}_R1_unpaired.fastq.gz" \
    "${PREFIX}_R2.fastq.gz" "${PREFIX}_R2_unpaired.fastq.gz" \
    SLIDINGWINDOW:4:20 \
    MINLEN:36

if [[ $? -eq 0 ]]; then
    echo -e "\n=== Output Files ==="
    ls -lh "${PREFIX}"*.fastq.gz

    echo -e "\n=== Read Counts ==="
    BEFORE=$(zcat "$R1" | wc -l | awk '{print $1/4}')
    AFTER=$(zcat "${PREFIX}_R1.fastq.gz" | wc -l | awk '{print $1/4}')
    UNPAIRED=$(zcat "${PREFIX}_R1_unpaired.fastq.gz" "${PREFIX}_R2_unpaired.fastq.gz" | wc -l | awk '{print $1/4}')

    echo "Input pairs: $BEFORE"
    echo "Output pairs: $AFTER"
    echo "Unpaired: $UNPAIRED"
    echo "Retention: $(echo "scale=1; $AFTER * 100 / $BEFORE" | bc)%"
else
    echo "Quality filtering failed"
    exit 1
fi

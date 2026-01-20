#!/bin/bash
# Process Hi-C pairs with pairtools

SAMPLE=${1:-sample}
CHROMSIZES=${2:-chromsizes.txt}

echo "Processing $SAMPLE..."

# Parse, sort, dedup, filter
pairtools parse \
    --chroms-path $CHROMSIZES \
    --min-mapq 30 \
    --walks-policy 5unique \
    ${SAMPLE}.bam | \
pairtools sort --nproc 4 | \
pairtools dedup --max-mismatch 1 --output-stats ${SAMPLE}.dedup_stats.txt | \
pairtools select '(pair_type == "UU")' \
    --output ${SAMPLE}.valid.pairs.gz

# Generate stats
pairtools stats ${SAMPLE}.valid.pairs.gz > ${SAMPLE}.stats.txt

echo "Stats:"
grep -E "^total|^cis|^trans" ${SAMPLE}.stats.txt

# Create cooler matrix
cooler cload pairs \
    -c1 2 -p1 3 -c2 4 -p2 5 \
    ${CHROMSIZES}:10000 \
    ${SAMPLE}.valid.pairs.gz \
    ${SAMPLE}.cool

echo "Done. Output: ${SAMPLE}.valid.pairs.gz, ${SAMPLE}.cool"

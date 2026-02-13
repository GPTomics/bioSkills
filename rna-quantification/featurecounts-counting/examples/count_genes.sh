#!/bin/bash
# Reference: DESeq2 1.42+, HISAT2 2.2.1+, STAR 2.7.11+, Subread 2.0+, edgeR 4.0+, pandas 2.2+, scanpy 1.10+ | Verify API if version differs
# Count genes from BAM files using featureCounts

GTF="Homo_sapiens.GRCh38.110.gtf"
OUTPUT="gene_counts.txt"
THREADS=8

# Paired-end, reverse stranded (Illumina TruSeq)
featureCounts \
    -p --countReadPairs \
    -s 2 \
    -T $THREADS \
    -a $GTF \
    -o $OUTPUT \
    *.bam

echo "Counts written to $OUTPUT"
echo "Summary statistics in ${OUTPUT}.summary"

# Extract clean count matrix
cut -f1,7- $OUTPUT | tail -n +2 > count_matrix.txt
echo "Clean matrix written to count_matrix.txt"

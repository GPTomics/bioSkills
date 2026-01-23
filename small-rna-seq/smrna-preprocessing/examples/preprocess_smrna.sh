#!/bin/bash
# Preprocess small RNA-seq data for miRNA analysis

INPUT=$1
OUTPUT_DIR=${2:-"preprocessed"}
ADAPTER=${3:-"TGGAATTCTCGGGTGCCAAGG"}  # Illumina TruSeq small RNA adapter

mkdir -p $OUTPUT_DIR

SAMPLE=$(basename $INPUT .fastq.gz)

echo "Processing: $SAMPLE"
echo "Adapter: $ADAPTER"

# Step 1: Trim adapters and filter by size
# -m 18: Minimum length - miRNAs are 18-25 nt, piRNAs 26-32 nt
# -M 35: Maximum length - excludes degradation products and longer ncRNAs
# --discard-untrimmed: Reads without adapter are likely not small RNA
# -q 20: Trim bases with quality < 20 from 3' end
cutadapt \
    -a $ADAPTER \
    -m 18 \
    -M 35 \
    -q 20 \
    --discard-untrimmed \
    -o ${OUTPUT_DIR}/${SAMPLE}_trimmed.fastq.gz \
    $INPUT \
    > ${OUTPUT_DIR}/${SAMPLE}_cutadapt.log 2>&1

echo "Trimming complete. Stats:"
grep -E "Total reads|Reads with adapters|Reads written" ${OUTPUT_DIR}/${SAMPLE}_cutadapt.log

# Step 2: Size selection for specific small RNA classes
# miRNA: 18-26 nt (peak at 21-23)
cutadapt \
    -m 18 -M 26 \
    -o ${OUTPUT_DIR}/${SAMPLE}_miRNA_size.fastq.gz \
    ${OUTPUT_DIR}/${SAMPLE}_trimmed.fastq.gz \
    > ${OUTPUT_DIR}/${SAMPLE}_size_select.log 2>&1

echo "Size selection complete"

# Step 3: Collapse identical reads (speeds up mapping)
# Output FASTA with counts in header
zcat ${OUTPUT_DIR}/${SAMPLE}_miRNA_size.fastq.gz | \
    awk 'NR%4==2' | \
    sort | uniq -c | sort -rn | \
    awk '{print ">"NR"_x"$1"\n"$2}' \
    > ${OUTPUT_DIR}/${SAMPLE}_collapsed.fa

UNIQUE=$(wc -l < ${OUTPUT_DIR}/${SAMPLE}_collapsed.fa)
UNIQUE=$((UNIQUE / 2))
echo "Unique sequences: $UNIQUE"

# Step 4: Length distribution for QC
echo "Length distribution:"
zcat ${OUTPUT_DIR}/${SAMPLE}_miRNA_size.fastq.gz | \
    awk 'NR%4==2 {lengths[length($0)]++} END {for (l in lengths) print l, lengths[l]}' | \
    sort -n

echo "Output files in $OUTPUT_DIR:"
ls -lh $OUTPUT_DIR/${SAMPLE}*

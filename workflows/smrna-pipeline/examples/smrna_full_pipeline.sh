#!/bin/bash
# Reference: cutadapt 4.4+, miRDeep2 2.0.1.3+, bowtie 1.3+ | Verify API if version differs
# Small RNA-seq pipeline -- miRDeep2 NOVEL-DISCOVERY route (genome + bowtie1 + miRDeep2).
# For the common known-miRNA + isomiR route, use miRge3.0 (see SKILL.md primary path).
set -e

FASTQ=$1
# miRDeep2.pl runs after a `cd` into the output dir, so resolve genome/miRBase inputs to
# ABSOLUTE paths up front; relative positional args would not resolve after the cd.
GENOME_FA=$(realpath "$2")
GENOME_INDEX=$3  # bowtie index prefix (used before the cd, so relative is fine)
MATURE_FA=$(realpath "$4")     # mature miRNAs from miRBase
HAIRPIN_FA=$(realpath "$5")    # hairpin sequences from miRBase
OUTPUT_DIR=${6:-"smrna_results"}
ADAPTER=${7:-"TGGAATTCTCGGGTGCCAAGG"}  # Illumina TruSeq small RNA adapter

mkdir -p ${OUTPUT_DIR}/{trimmed,mirdeep2,counts}

echo "=== Step 1: Adapter Trimming ==="
# Remove adapter and filter by size
# 18-30nt: typical range for miRNAs (peak at 21-23nt)
cutadapt \
    -a $ADAPTER \
    --minimum-length 18 \
    --maximum-length 30 \
    --discard-untrimmed \
    -o ${OUTPUT_DIR}/trimmed/trimmed.fastq.gz \
    $FASTQ \
    > ${OUTPUT_DIR}/trimmed/cutadapt_report.txt

# QC: Check size distribution
echo "Size distribution after trimming:"
zcat ${OUTPUT_DIR}/trimmed/trimmed.fastq.gz | \
    awk 'NR%4==2 {print length}' | sort -n | uniq -c | head -20

echo "=== Step 2: miRDeep2 Alignment ==="
# Collapse identical reads and align to genome
# -e: FASTQ input
# -h: parse to FASTA
# -i: convert RNA to DNA
# -j: remove reads with non-canonical letters
# -l 18: minimum read length
# -m: collapse reads
# mapper.pl reads plain-text FASTQ (no gzip handling); decompress first.
gunzip -kc ${OUTPUT_DIR}/trimmed/trimmed.fastq.gz > ${OUTPUT_DIR}/trimmed/trimmed.fastq
mapper.pl ${OUTPUT_DIR}/trimmed/trimmed.fastq \
    -e -h -i -j -l 18 -m \
    -p $GENOME_INDEX \
    -s ${OUTPUT_DIR}/mirdeep2/reads_collapsed.fa \
    -t ${OUTPUT_DIR}/mirdeep2/reads_vs_genome.arf \
    -v

echo "=== Step 3: miRDeep2 Quantification ==="
cd "${OUTPUT_DIR}/mirdeep2"

# Run miRDeep2
# Outputs: known miRNA counts, novel miRNA predictions
miRDeep2.pl \
    reads_collapsed.fa \
    $GENOME_FA \
    reads_vs_genome.arf \
    $MATURE_FA \
    none \
    $HAIRPIN_FA \
    -t Human \
    2> mirdeep2_run.log

cd -

echo "=== Step 4: Extract Counts ==="
# Parse miRDeep2 output for count matrix
# Use miRNAs_expressed_all_samples file
cp ${OUTPUT_DIR}/mirdeep2/miRNAs_expressed_all_samples*.csv \
   ${OUTPUT_DIR}/counts/mirna_counts.csv

echo "=== Pipeline Complete ==="
echo "Results in: $OUTPUT_DIR"
echo "Key outputs:"
echo "  - Trimmed reads: ${OUTPUT_DIR}/trimmed/"
echo "  - miRNA counts: ${OUTPUT_DIR}/counts/mirna_counts.csv"
echo "  - Novel miRNAs: ${OUTPUT_DIR}/mirdeep2/result*.html"
echo ""
echo "Next steps:"
echo "  - Run DESeq2 on counts for differential expression"
echo "  - Run target prediction on significant miRNAs"

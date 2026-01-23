#!/bin/bash
# Ribo-seq preprocessing pipeline

INPUT=$1
OUTPUT_PREFIX=${2:-"riboseq"}
ADAPTER=${3:-"CTGTAGGCACCATCAAT"}  # Common Ribo-seq adapter

# Directories
mkdir -p logs

echo "Ribo-seq preprocessing: $INPUT"

# Step 1: Trim adapters
echo "Trimming adapters..."
cutadapt \
    -a $ADAPTER \
    -m 20 \
    -M 40 \
    -j 4 \
    -o ${OUTPUT_PREFIX}_trimmed.fastq.gz \
    $INPUT \
    > logs/${OUTPUT_PREFIX}_cutadapt.log 2>&1

# Step 2: Size selection
# Ribosome footprints are typically 28-32 nt
# Shorter reads may be degradation; longer may include linker
MIN_LEN=28  # Protected fragment minimum
MAX_LEN=32  # Protected fragment maximum
echo "Selecting ${MIN_LEN}-${MAX_LEN} nt reads..."
cutadapt \
    -m $MIN_LEN \
    -M $MAX_LEN \
    -o ${OUTPUT_PREFIX}_sized.fastq.gz \
    ${OUTPUT_PREFIX}_trimmed.fastq.gz \
    > logs/${OUTPUT_PREFIX}_size_select.log 2>&1

SIZED_READS=$(zcat ${OUTPUT_PREFIX}_sized.fastq.gz | wc -l)
SIZED_READS=$((SIZED_READS / 4))
echo "Reads in size range: $SIZED_READS"

# Step 3: rRNA removal
# Critical step - Ribo-seq can have 50-90% rRNA contamination
echo "Removing rRNA..."
if command -v sortmerna &> /dev/null; then
    # SortMeRNA (more comprehensive)
    sortmerna \
        --ref /path/to/rRNA_databases/silva-euk-18s-id95.fasta \
        --ref /path/to/rRNA_databases/silva-euk-28s-id98.fasta \
        --reads ${OUTPUT_PREFIX}_sized.fastq.gz \
        --aligned ${OUTPUT_PREFIX}_rRNA \
        --other ${OUTPUT_PREFIX}_non_rRNA \
        --fastx \
        --threads 8 \
        > logs/${OUTPUT_PREFIX}_sortmerna.log 2>&1
    NON_RRNA="${OUTPUT_PREFIX}_non_rRNA.fastq.gz"
else
    # Bowtie2 fallback
    bowtie2 -x rRNA_index \
        -U ${OUTPUT_PREFIX}_sized.fastq.gz \
        --un ${OUTPUT_PREFIX}_non_rRNA.fastq.gz \
        -S /dev/null \
        -p 8 \
        2> logs/${OUTPUT_PREFIX}_rrna_removal.log
    NON_RRNA="${OUTPUT_PREFIX}_non_rRNA.fastq.gz"
fi

# Step 4: Align to transcriptome
echo "Aligning to transcriptome..."
STAR --runMode alignReads \
    --genomeDir /path/to/STAR_index \
    --readFilesIn $NON_RRNA \
    --readFilesCommand zcat \
    --outFilterMultimapNmax 1 \
    --outFilterMismatchNmax 2 \
    --alignIntronMax 1 \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix ${OUTPUT_PREFIX}_ \
    --runThreadN 8 \
    > logs/${OUTPUT_PREFIX}_star.log 2>&1

# Index BAM
samtools index ${OUTPUT_PREFIX}_Aligned.sortedByCoord.out.bam

# Summary statistics
echo ""
echo "Preprocessing complete!"
echo "Read length distribution:"
samtools view ${OUTPUT_PREFIX}_Aligned.sortedByCoord.out.bam | \
    awk '{lengths[length($10)]++} END {for (l in lengths) print l, lengths[l]}' | \
    sort -k1n | head -10

echo ""
echo "Alignment stats:"
samtools flagstat ${OUTPUT_PREFIX}_Aligned.sortedByCoord.out.bam

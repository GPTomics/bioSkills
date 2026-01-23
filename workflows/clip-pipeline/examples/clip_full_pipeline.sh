#!/bin/bash
# Complete CLIP-seq analysis pipeline

FASTQ=$1
STAR_INDEX=$2
GENOME_FA=$3
GTF=$4
OUTPUT_DIR=${5:-"clip_results"}
UMI_PATTERN=${6:-"NNNNNNNNNN"}  # 10nt UMI at read start
ADAPTER=${7:-"AGATCGGAAGAGCACACGTCT"}
SPECIES=${8:-"hg38"}
THREADS=${9:-8}

mkdir -p ${OUTPUT_DIR}/{preprocessed,aligned,peaks,motifs}

echo "=== Step 1: UMI Extraction ==="
# Extract UMI from read start
umi_tools extract \
    --stdin=${FASTQ} \
    --bc-pattern=${UMI_PATTERN} \
    --stdout=${OUTPUT_DIR}/preprocessed/umi_extracted.fastq.gz \
    --log=${OUTPUT_DIR}/preprocessed/umi_extract.log

echo "=== Step 2: Adapter Trimming ==="
cutadapt \
    -a ${ADAPTER} \
    --minimum-length 20 \
    --quality-cutoff 20 \
    -o ${OUTPUT_DIR}/preprocessed/trimmed.fastq.gz \
    ${OUTPUT_DIR}/preprocessed/umi_extracted.fastq.gz \
    > ${OUTPUT_DIR}/preprocessed/cutadapt_report.txt

echo "=== Step 3: Alignment ==="
STAR \
    --genomeDir ${STAR_INDEX} \
    --readFilesIn ${OUTPUT_DIR}/preprocessed/trimmed.fastq.gz \
    --readFilesCommand zcat \
    --outFilterMismatchNmax 2 \
    --outFilterMultimapNmax 1 \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix ${OUTPUT_DIR}/aligned/ \
    --runThreadN ${THREADS}

samtools index ${OUTPUT_DIR}/aligned/Aligned.sortedByCoord.out.bam

echo "=== Step 4: UMI Deduplication ==="
umi_tools dedup \
    -I ${OUTPUT_DIR}/aligned/Aligned.sortedByCoord.out.bam \
    -S ${OUTPUT_DIR}/aligned/dedup.bam \
    --output-stats=${OUTPUT_DIR}/aligned/dedup_stats \
    --log=${OUTPUT_DIR}/aligned/dedup.log

samtools index ${OUTPUT_DIR}/aligned/dedup.bam

# QC: Check duplication rate
echo "Deduplication stats:"
grep "Input Reads" ${OUTPUT_DIR}/aligned/dedup.log
grep "Output Reads" ${OUTPUT_DIR}/aligned/dedup.log

echo "=== Step 5: Peak Calling ==="
clipper \
    -b ${OUTPUT_DIR}/aligned/dedup.bam \
    -s ${SPECIES} \
    -o ${OUTPUT_DIR}/peaks/peaks.bed \
    --save-pickle

echo "Peaks called: $(wc -l < ${OUTPUT_DIR}/peaks/peaks.bed)"

echo "=== Step 6: Peak Annotation ==="
# Annotate to genomic features
bedtools intersect \
    -a ${OUTPUT_DIR}/peaks/peaks.bed \
    -b ${GTF} \
    -wo \
    > ${OUTPUT_DIR}/peaks/peaks_annotated.txt

echo "=== Step 7: Motif Analysis ==="
# Extract peak sequences
bedtools getfasta \
    -fi ${GENOME_FA} \
    -bed ${OUTPUT_DIR}/peaks/peaks.bed \
    -fo ${OUTPUT_DIR}/motifs/peaks.fa

# HOMER de novo motif discovery
# -rna: RNA mode (U instead of T)
# -len 6,7,8: typical RBP motif lengths
findMotifs.pl \
    ${OUTPUT_DIR}/motifs/peaks.fa \
    fasta \
    ${OUTPUT_DIR}/motifs/ \
    -rna \
    -len 6,7,8

echo "=== Pipeline Complete ==="
echo "Results in: ${OUTPUT_DIR}"
echo ""
echo "Key outputs:"
echo "  - Deduplicated BAM: ${OUTPUT_DIR}/aligned/dedup.bam"
echo "  - Binding sites: ${OUTPUT_DIR}/peaks/peaks.bed"
echo "  - Motifs: ${OUTPUT_DIR}/motifs/homerResults.html"

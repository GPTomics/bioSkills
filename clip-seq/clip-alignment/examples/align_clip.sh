#!/bin/bash
# CLIP-seq alignment

FASTQ=$1
STAR_INDEX=$2
OUTPUT_PREFIX=${3:-"clip"}

STAR --runMode alignReads \
    --genomeDir $STAR_INDEX \
    --readFilesIn $FASTQ \
    --readFilesCommand zcat \
    --outFilterMultimapNmax 1 \
    --outFilterMismatchNmax 1 \
    --alignEndsType EndToEnd \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix ${OUTPUT_PREFIX}_ \
    --runThreadN 8

samtools index ${OUTPUT_PREFIX}_Aligned.sortedByCoord.out.bam

# Deduplicate with UMIs
umi_tools dedup \
    --stdin=${OUTPUT_PREFIX}_Aligned.sortedByCoord.out.bam \
    --stdout=${OUTPUT_PREFIX}_deduped.bam

samtools index ${OUTPUT_PREFIX}_deduped.bam

echo "Alignment complete: ${OUTPUT_PREFIX}_deduped.bam"

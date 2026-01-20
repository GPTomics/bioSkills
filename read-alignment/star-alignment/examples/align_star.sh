#!/bin/bash
# RNA-seq alignment with STAR

set -euo pipefail

R1=${1:-reads_1.fastq.gz}
R2=${2:-reads_2.fastq.gz}
GENOME_DIR=${3:-star_index}
OUTPUT_PREFIX=${4:-star_out/sample}
THREADS=${5:-8}

mkdir -p $(dirname $OUTPUT_PREFIX)

echo "=== STAR RNA-seq Alignment ==="

STAR --runThreadN $THREADS \
    --genomeDir $GENOME_DIR \
    --readFilesIn $R1 $R2 \
    --readFilesCommand zcat \
    --outFileNamePrefix ${OUTPUT_PREFIX}_ \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --quantMode GeneCounts \
    --twopassMode Basic

# Index BAM
samtools index ${OUTPUT_PREFIX}_Aligned.sortedByCoord.out.bam

echo "BAM: ${OUTPUT_PREFIX}_Aligned.sortedByCoord.out.bam"
echo "Gene counts: ${OUTPUT_PREFIX}_ReadsPerGene.out.tab"

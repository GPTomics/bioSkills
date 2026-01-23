#!/bin/bash
# MeRIP-seq alignment pipeline

GENOME_DIR=$1
SAMPLE_SHEET=$2
OUTPUT_DIR=${3:-"aligned"}
THREADS=${4:-8}

mkdir -p $OUTPUT_DIR

# Sample sheet format: sample_name,read1,read2,sample_type (IP or Input)
while IFS=',' read -r sample r1 r2 sample_type; do
    echo "Aligning $sample ($sample_type)..."

    STAR --genomeDir $GENOME_DIR \
        --readFilesIn $r1 $r2 \
        --readFilesCommand zcat \
        --runThreadN $THREADS \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${OUTPUT_DIR}/${sample}_ \
        --outSAMattributes Standard

    # Index BAM
    samtools index ${OUTPUT_DIR}/${sample}_Aligned.sortedByCoord.out.bam

    # Basic stats
    samtools flagstat ${OUTPUT_DIR}/${sample}_Aligned.sortedByCoord.out.bam \
        > ${OUTPUT_DIR}/${sample}_flagstat.txt

done < $SAMPLE_SHEET

echo "Alignment complete. BAMs in $OUTPUT_DIR"

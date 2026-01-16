#!/bin/bash

GENOME_DIR="genome"
READS_DIR="fastq"
OUTPUT_DIR="aligned"

mkdir -p $OUTPUT_DIR

bismark_genome_preparation --bowtie2 $GENOME_DIR

bismark --genome $GENOME_DIR \
    -1 ${READS_DIR}/sample_R1.fastq.gz \
    -2 ${READS_DIR}/sample_R2.fastq.gz \
    --parallel 4 \
    -o $OUTPUT_DIR

deduplicate_bismark --paired \
    --bam ${OUTPUT_DIR}/sample_R1_bismark_bt2_pe.bam

samtools sort ${OUTPUT_DIR}/sample_R1_bismark_bt2_pe.deduplicated.bam \
    -o ${OUTPUT_DIR}/sample.sorted.bam
samtools index ${OUTPUT_DIR}/sample.sorted.bam

cat ${OUTPUT_DIR}/*_report.txt

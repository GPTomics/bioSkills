#!/bin/bash

GENOME_DIR="genome"
BAM_DIR="aligned"
OUTPUT_DIR="methylation"

mkdir -p $OUTPUT_DIR

for bam in ${BAM_DIR}/*.deduplicated.bam; do
    sample=$(basename $bam .deduplicated.bam)

    bismark_methylation_extractor \
        --paired-end \
        --no_overlap \
        --gzip \
        --bedGraph \
        --cytosine_report \
        --genome_folder $GENOME_DIR \
        --parallel 4 \
        -o $OUTPUT_DIR \
        $bam
done

ls -la ${OUTPUT_DIR}/*.bismark.cov.gz

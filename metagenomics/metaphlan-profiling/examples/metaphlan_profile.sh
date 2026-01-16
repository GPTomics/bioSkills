#!/bin/bash

READS_DIR="fastq"
OUTPUT_DIR="metaphlan_output"
NPROC=8

mkdir -p $OUTPUT_DIR/profiles $OUTPUT_DIR/mapout

for fq in ${READS_DIR}/*.fastq.gz; do
    sample=$(basename $fq .fastq.gz)
    echo "Processing ${sample}..."

    metaphlan $fq \
        --input_type fastq \
        --nproc $NPROC \
        --output_file ${OUTPUT_DIR}/profiles/${sample}_profile.txt \
        --mapout ${OUTPUT_DIR}/mapout/${sample}.map.bz2
done

echo "Merging profiles..."
merge_metaphlan_tables.py ${OUTPUT_DIR}/profiles/*_profile.txt > ${OUTPUT_DIR}/merged_abundance.txt

echo "Top species across all samples:"
grep "s__" ${OUTPUT_DIR}/merged_abundance.txt | grep -v "t__" | head -20

echo ""
echo "Results saved to ${OUTPUT_DIR}/merged_abundance.txt"

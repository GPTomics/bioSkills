#!/bin/bash

READS="reads.fastq.gz"
OUTPUT_DIR="qc_output"
THREADS=4

mkdir -p $OUTPUT_DIR

echo "Generating basic statistics..."
NanoStat --fastq $READS --threads $THREADS > ${OUTPUT_DIR}/nanostats.txt
cat ${OUTPUT_DIR}/nanostats.txt

echo ""
echo "Generating QC plots..."
NanoPlot --fastq $READS \
    -o ${OUTPUT_DIR}/nanoplot \
    -t $THREADS \
    --N50 \
    --title "Read QC" \
    --plots hex dot \
    --format png

echo "QC report saved to ${OUTPUT_DIR}/nanoplot/NanoPlot-report.html"

MIN_QUALITY=10
MIN_LENGTH=1000

echo ""
echo "Filtering reads (Q>=${MIN_QUALITY}, L>=${MIN_LENGTH})..."
gunzip -c $READS | chopper \
    --quality $MIN_QUALITY \
    --minlength $MIN_LENGTH \
    --threads $THREADS \
    | gzip > ${OUTPUT_DIR}/filtered.fastq.gz

echo ""
echo "Filtered read statistics:"
NanoStat --fastq ${OUTPUT_DIR}/filtered.fastq.gz --threads $THREADS

echo ""
echo "Summary:"
original=$(seqkit stats $READS -T | tail -1 | cut -f 4)
filtered=$(seqkit stats ${OUTPUT_DIR}/filtered.fastq.gz -T | tail -1 | cut -f 4)
echo "Original reads: $original"
echo "Filtered reads: $filtered"
pct=$(echo "scale=1; $filtered * 100 / $original" | bc)
echo "Retained: ${pct}%"

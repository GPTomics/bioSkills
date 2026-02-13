#!/bin/bash
# Reference: MACS3 3.0+, Subread 2.0+, bedtools 2.31+, deepTools 3.5+, pybedtools 0.9+, pysam 0.22+, samtools 1.19+ | Verify API if version differs
# Run IDR analysis on two replicates

REP1=$1
REP2=$2
OUTPUT=$3

if [ -z "$REP1" ] || [ -z "$REP2" ]; then
    echo "Usage: $0 <rep1_peaks> <rep2_peaks> [output_prefix]"
    exit 1
fi

OUTPUT=${OUTPUT:-idr_output}

sort -k7,7nr "$REP1" > rep1_sorted.narrowPeak
sort -k7,7nr "$REP2" > rep2_sorted.narrowPeak

idr --samples rep1_sorted.narrowPeak rep2_sorted.narrowPeak \
    --input-file-type narrowPeak \
    --rank signal.value \
    --output-file "${OUTPUT}.txt" \
    --plot "${OUTPUT}.pdf" \
    --log-output-file "${OUTPUT}.log"

echo "Peaks at IDR < 0.05: $(awk '$5 >= 540' ${OUTPUT}.txt | wc -l)"
echo "Peaks at IDR < 0.1: $(awk '$5 >= 415' ${OUTPUT}.txt | wc -l)"

rm -f rep1_sorted.narrowPeak rep2_sorted.narrowPeak

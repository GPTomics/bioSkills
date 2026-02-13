#!/bin/bash
# Reference: MACS3 3.0+, bedtools 2.31+, pandas 2.2+, samtools 1.19+ | Verify API if version differs
# CLIP-seq peak calling

BAM=$1
SPECIES=${2:-"hg38"}
OUTPUT=${3:-"peaks.bed"}

clipper -b $BAM -s $SPECIES -o $OUTPUT --save-pickle

echo "Peaks called: $OUTPUT"
wc -l $OUTPUT

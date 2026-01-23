#!/bin/bash
# CLIP-seq peak calling

BAM=$1
SPECIES=${2:-"hg38"}
OUTPUT=${3:-"peaks.bed"}

clipper -b $BAM -s $SPECIES -o $OUTPUT --save-pickle

echo "Peaks called: $OUTPUT"
wc -l $OUTPUT

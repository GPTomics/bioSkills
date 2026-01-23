#!/bin/bash
# CLIP-seq motif analysis

PEAKS=$1
GENOME=$2
OUTPUT_DIR=${3:-"motif_output"}

mkdir -p $OUTPUT_DIR

# Extract peak sequences
bedtools getfasta -fi $GENOME -bed $PEAKS -fo ${OUTPUT_DIR}/peaks.fa

# HOMER de novo motif discovery
# -rna: RNA motifs (use U instead of T)
# -len: motif lengths to search
findMotifs.pl ${OUTPUT_DIR}/peaks.fa fasta $OUTPUT_DIR \
    -len 6,7,8 \
    -rna

echo "Motif analysis complete: $OUTPUT_DIR"

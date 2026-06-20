#!/bin/bash
# Reference: Ribo-TISH 0.2.7+ | Verify API if version differs
# Map translation initiation sites from initiation-drug (TI-seq) ribosome profiling.

set -euo pipefail

ELONG_BAM=$1     # standard elongation Ribo-seq BAM
TIS_BAM=$2       # initiation-drug BAM (harringtonine or LTM)
GTF=$3
GENOME=$4
OUTPREFIX=${5:-tis}
DRUG=${6:-harr}  # "harr" for harringtonine-type libraries; "" for LTM

# Step 1: QC each library. ribotish quality writes a <bam>.para.py offset file
# and a QC figure; confirms the drug enriched start-codon signal.
ribotish quality -b "$ELONG_BAM" -g "$GTF" -o "${OUTPREFIX}_elong_quality.txt"
ribotish quality -b "$TIS_BAM" -g "$GTF" -o "${OUTPREFIX}_tis_quality.txt"

# Step 2: predict initiation sites. --alt enables near-cognate (CUG/GUG/...) starts,
# which most uORFs use. --harr marks a harringtonine-type TIS library (broader peaks).
HARR_ARG=""
if [ "$DRUG" = "harr" ]; then
    HARR_ARG="--harr --harrwidth 15"
fi
ribotish predict \
    -b "$ELONG_BAM" \
    -t "$TIS_BAM" \
    -g "$GTF" \
    -f "$GENOME" \
    $HARR_ARG --alt \
    -o "${OUTPREFIX}_predictions.txt"

echo "Initiation sites written to ${OUTPREFIX}_predictions.txt"
echo "Columns include start codon, ORF type, and significance."
echo "For cryptic/near-cognate EM detection, also consider PRICE:"
echo "  gedi -e Price -reads $ELONG_BAM -genomic <prepared_genome> -prefix ${OUTPREFIX}_price"

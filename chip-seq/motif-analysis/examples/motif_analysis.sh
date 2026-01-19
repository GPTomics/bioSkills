#!/bin/bash
# ChIP-seq motif analysis with HOMER and MEME
set -euo pipefail

PEAKS=$1
GENOME=$2  # hg38 or mm10
OUTDIR=$3
GENOME_FA=${4:-${GENOME}.fa}

mkdir -p ${OUTDIR}/homer ${OUTDIR}/meme

echo "=== Motif Analysis Pipeline ==="
echo "Peaks: $PEAKS"
echo "Genome: $GENOME"

# Count peaks
NPEAK=$(wc -l < $PEAKS)
echo "Number of peaks: $NPEAK"

# HOMER de novo and known motif analysis
echo "Running HOMER..."
findMotifsGenome.pl $PEAKS $GENOME ${OUTDIR}/homer \
    -size 200 \
    -mask \
    -p 4

# Extract sequences for MEME (center +/- 100bp)
echo "Extracting peak sequences..."
awk 'BEGIN{OFS="\t"} {
    center = int(($2+$3)/2)
    start = center - 100
    if(start < 0) start = 0
    print $1, start, center + 100
}' $PEAKS > ${OUTDIR}/peaks_centered.bed

bedtools getfasta \
    -fi $GENOME_FA \
    -bed ${OUTDIR}/peaks_centered.bed \
    -fo ${OUTDIR}/peaks.fa

# MEME de novo discovery
echo "Running MEME..."
meme ${OUTDIR}/peaks.fa \
    -dna \
    -oc ${OUTDIR}/meme \
    -mod zoops \
    -nmotifs 10 \
    -minw 6 \
    -maxw 15 \
    -revcomp

# Summary
echo ""
echo "=== Results ==="
echo "HOMER results: ${OUTDIR}/homer/homerResults.html"
echo "HOMER known: ${OUTDIR}/homer/knownResults.txt"
echo "MEME results: ${OUTDIR}/meme/meme.html"

# Print top HOMER known motifs
echo ""
echo "Top 10 known motifs (HOMER):"
head -11 ${OUTDIR}/homer/knownResults.txt | tail -10 | cut -f1,3,6,7

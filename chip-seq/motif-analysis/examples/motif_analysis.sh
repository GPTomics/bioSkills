#!/bin/bash
# Reference: HOMER 4.11+, MEME suite 5.5+ (STREME replaces DREME), JASPAR 2024 CORE, bedtools 2.31+ | Verify API if version differs
# ChIP-seq motif analysis: HOMER de novo + known, MEME-ChIP comprehensive
# (STREME + CentriMo + TOMTOM + FIMO), AME differential. Uses summit-centered
# ±100 bp sequences (motif enrichment improves dramatically vs full peak width).

set -euo pipefail

NARROWPEAK=$1   # MACS narrowPeak format (column 10 = summit offset from start)
GENOME=$2       # hg38 / mm10
OUTDIR=$3
GENOME_FA=${4:-${GENOME}.fa}
JASPAR_MEME=${5:-JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt}

mkdir -p ${OUTDIR}/{homer,meme_chip,fimo}

NPEAK=$(wc -l < $NARROWPEAK)
echo "=== Motif Analysis: $NPEAK peaks ==="
[ $NPEAK -lt 500 ] && echo "WARNING: <500 peaks; motif discovery may be underpowered"


# === Recenter on summit ±100 bp (narrowPeak col 10 is summit offset) ===
# Motif enrichment improves substantially when sequences are summit-centered
# (200 bp window) vs full peak width. CentriMo requires this for central
# enrichment testing.
awk 'BEGIN{OFS="\t"} {
    summit = $2 + $10
    start = summit - 100
    if (start < 0) start = 0
    end = summit + 100
    print $1, start, end, $4, $5, $6
}' $NARROWPEAK > ${OUTDIR}/peaks_summit_centered.bed

bedtools getfasta -fi $GENOME_FA -bed ${OUTDIR}/peaks_summit_centered.bed \
    -fo ${OUTDIR}/peaks_summit_centered.fa


# === HOMER: de novo + known motif discovery ===
# -mask = mask repeats; prevents AluY/LINE/LTR consensus from dominating de novo
# -p 8 = parallel; -size 200 = use ±100 bp window we already prepared
# Background: HOMER auto-samples GC-matched genomic regions. For large peak sets
# (>5% of genome) supply -bg with explicit random/control intervals.
findMotifsGenome.pl ${OUTDIR}/peaks_summit_centered.bed $GENOME ${OUTDIR}/homer \
    -size 200 \
    -mask \
    -p 8 \
    -len 8,10,12


# === MEME-ChIP: comprehensive pipeline ===
# Runs STREME (de novo, replaces DREME from 5.4+), MEME (longer / gapped),
# CentriMo (central enrichment of known motifs — direct vs tethered binding),
# TOMTOM (compare to JASPAR), FIMO (motif instance scanning).
# Markov order-2 background is the internal default; preserves trinucleotide
# composition (CpG context).
meme-chip \
    -oc ${OUTDIR}/meme_chip \
    -db $JASPAR_MEME \
    -meme-nmotifs 5 -meme-minw 6 -meme-maxw 20 \
    -streme-nmotifs 10 -streme-minw 4 -streme-maxw 15 \
    -centrimo-score 5.0 -centrimo-ethresh 10.0 \
    ${OUTDIR}/peaks_summit_centered.fa


# === FIMO: genome-wide motif scanning at p ≤ 1e-5 ===
# Default p ≤ 1e-4 produces millions of false positives genome-wide; tighten
# for whole-genome OR restrict to peaks for finer p-value.
fimo --oc ${OUTDIR}/fimo \
    --thresh 1e-5 \
    --max-stored-scores 1000000 \
    $JASPAR_MEME \
    ${OUTDIR}/peaks_summit_centered.fa


# === Optional: AME differential between two peak sets ===
# Use when comparing peaks gained in condition A vs gained in condition B.
# Requires peaks_A.fa and peaks_B.fa; uncomment + supply paths.
# ame --oc ${OUTDIR}/ame_diff \
#     --control peaks_B.fa \
#     --scoring avg --method fisher \
#     peaks_A.fa $JASPAR_MEME


# === Summary ===
echo ""
echo "=== Results ==="
echo "HOMER de novo + known: ${OUTDIR}/homer/homerResults.html"
echo "HOMER known table: ${OUTDIR}/homer/knownResults.txt"
echo "MEME-ChIP report: ${OUTDIR}/meme_chip/meme-chip.html"
echo "FIMO instances: ${OUTDIR}/fimo/fimo.tsv"

# Print top 10 known motifs from HOMER
echo ""
echo "Top 10 HOMER known motifs:"
if [ -f ${OUTDIR}/homer/knownResults.txt ]; then
    head -11 ${OUTDIR}/homer/knownResults.txt | tail -10 | \
        awk -F'\t' '{printf "%-30s\tp=%s\ttarget=%s\tbg=%s\n", $1, $3, $6, $7}'
fi

# Sanity check: known TF should appear in top results for validated ChIP.
# Failure to recover known motif is a QC failure indicator — re-check FRiP,
# antibody validation, hyper-ChIPable contamination.

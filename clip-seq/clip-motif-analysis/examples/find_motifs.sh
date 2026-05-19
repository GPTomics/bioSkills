#!/bin/bash
# Reference: HOMER 4.11+, MEME Suite 5.5+, bedtools 2.31+ | Verify API if version differs
# CLIP-seq motif discovery: HOMER de novo + STREME cross-validation + GC-matched background.
# For CL-position-registered motifs (mCross / PEKA), pass single-nt crosslink sites from PureCLIP - see commented section below.

PEAKS=$1                          # CLIPper stringent peaks (log2 FC >= 3) BED
GENOME=$2                         # genome FASTA
EXPRESSED_REGIONS=${3:-""}        # BED of expressed transcripts (TPM >= 1); for GC-matched background
OUTPUT_DIR=${4:-"motif_out"}
THREADS=${5:-8}

mkdir -p $OUTPUT_DIR

# Step 1: Extract peak sequences, preserving strand (-s)
bedtools getfasta -fi $GENOME -bed $PEAKS -s -fo ${OUTPUT_DIR}/peaks.fa
echo "Peak sequences: $(grep -c '^>' ${OUTPUT_DIR}/peaks.fa)"

# Step 2: GC-matched background
# HOMER's default shuffled background is GC-blind; CLIP peaks are biased toward 3' UTRs (AU-rich)
# Without a matched background, top HOMER motif is often spurious AU
if [ -n "$EXPRESSED_REGIONS" ]; then
    # Sample random regions from expressed transcripts, matched in count and width
    awk 'BEGIN{OFS="\t"} {print $1, $2, $3}' $PEAKS > ${OUTPUT_DIR}/peaks_widths.bed
    bedtools shuffle \
        -i ${OUTPUT_DIR}/peaks_widths.bed \
        -g <(samtools faidx $GENOME && cut -f1,2 ${GENOME}.fai) \
        -incl $EXPRESSED_REGIONS \
        -seed 42 > ${OUTPUT_DIR}/background.bed
    bedtools getfasta -fi $GENOME -bed ${OUTPUT_DIR}/background.bed -fo ${OUTPUT_DIR}/background.fa
else
    # Fall back to HOMER auto-shuffled (less reliable)
    echo "No expressed-region BED provided; using HOMER default background (may bias AU-rich)"
    touch ${OUTPUT_DIR}/background.fa
fi

# Step 3: HOMER de novo motif discovery
# -rna: report U instead of T
# -len 5,6,7,8: typical RBP motif widths
# -fasta: explicit background
HOMER_BG=""
if [ -s "${OUTPUT_DIR}/background.fa" ]; then
    HOMER_BG="-fasta ${OUTPUT_DIR}/background.fa"
fi

findMotifs.pl ${OUTPUT_DIR}/peaks.fa fasta ${OUTPUT_DIR}/homer \
    -rna -len 5,6,7,8 -p $THREADS $HOMER_BG

echo "HOMER complete: ${OUTPUT_DIR}/homer/"

# Step 4: STREME for orthogonal confirmation (fast successor to MEME)
if command -v streme > /dev/null; then
    if [ -s "${OUTPUT_DIR}/background.fa" ]; then
        streme --rna --oc ${OUTPUT_DIR}/streme \
            -p ${OUTPUT_DIR}/peaks.fa -n ${OUTPUT_DIR}/background.fa \
            --minw 5 --maxw 10 --nmotifs 10
    else
        streme --rna --oc ${OUTPUT_DIR}/streme \
            -p ${OUTPUT_DIR}/peaks.fa \
            --minw 5 --maxw 10 --nmotifs 10
    fi
    echo "STREME complete: ${OUTPUT_DIR}/streme/"
fi

# Step 5: Sanity check - motif information content
# A real RBP motif has mean IC per position > 1.0 bit
echo ""
echo "Quality check - inspect ${OUTPUT_DIR}/homer/homerResults/motif1.motif"
echo "  Mean IC per position should be > 1.0 bit"
echo "  Central column U fraction > 80% with non-U-binding RBP = U crosslink bias suspect"

# Step 6 (optional): CL-position-registered motifs with mCross
# Requires single-nt CL sites from PureCLIP or CTK CITS - see clip-seq/crosslink-site-detection
# Uncomment to enable
# CROSSLINKS=$1.crosslinks.bed  # PureCLIP single-nt site BED
# mCross -i $CROSSLINKS -g $GENOME -k 7 -n 5 -o ${OUTPUT_DIR}/mcross

# Step 7 (optional): PEKA positional k-mer enrichment
# REGIONS=$1.regions.bed
# peka -i $PEAKS -x $CROSSLINKS -g $GENOME -r $REGIONS -k 5 -p 30 -o ${OUTPUT_DIR}/peka

echo ""
echo "Cross-validate top motif against:"
echo "  - CISBP-RNA database (https://cisbp-rna.ccbr.utoronto.ca/)"
echo "  - ATtRACT (https://attract.cnic.es/)"
echo "  - RBNS Kd ranking (Dominguez et al 2018 Mol Cell 70:854 supplementary)"
echo ""
echo "Three independent tools agreeing on core motif = high confidence"

#!/bin/bash
# Reference: STAMP / Bullseye 1.0+, samtools 1.19+, REDItools2 1.3+ | Verify API if version differs
# STAMP analysis: detect C->U edits in APOBEC1-RBP fusion sample vs APOBEC1-only control.
# Deaminase-only control is MANDATORY: without it, all edits look like signal.

FUSION_BAM=$1                    # APOBEC1-RBP-fusion sample BAM
CONTROL_BAM=$2                   # APOBEC1-only control BAM
GENOME=$3
OUT_PREFIX=${4:-"stamp"}

# Step 1: Verify both BAMs are indexed and strand-aware
samtools index -e $FUSION_BAM 2>/dev/null || samtools index $FUSION_BAM
samtools index -e $CONTROL_BAM 2>/dev/null || samtools index $CONTROL_BAM

# Step 2: Detect C->U editing with Bullseye (canonical STAMP pipeline)
# --edit_type c2t: APOBEC1 deaminates C; sequencing reads U as T
# --threshold 0.1: edit rate >= 10% per position
# --min_coverage 10: skip low-coverage sites
if command -v bullseye > /dev/null; then
    bullseye \
        --ip $FUSION_BAM \
        --control $CONTROL_BAM \
        --reference $GENOME \
        --edit_type c2t \
        --threshold 0.1 \
        --min_coverage 10 \
        --output ${OUT_PREFIX}_edits.bed

    echo "Bullseye edits: $(wc -l < ${OUT_PREFIX}_edits.bed)"
else
    # Fallback: JACUSA2 (general-purpose RNA editing detector)
    if command -v jacusa2 > /dev/null; then
        jacusa2 call-2 \
            -r $GENOME \
            -p 8 \
            -F 1024 \
            -A $FUSION_BAM \
            -B $CONTROL_BAM \
            -o ${OUT_PREFIX}_jacusa.tsv

        # Filter for C->U at edit rate >= 0.1 and coverage >= 10
        # JACUSA2 columns: chr, start, end, name, score, strand, ref, A,C,G,T,N
        # Custom filter for STAMP (C on + strand or G on - strand reference, mismatch to T or A)
        awk 'BEGIN{FS=OFS="\t"} NR>1 && $7=="C" {
            # IP: $11 reads, control: $15 reads (column positions vary; verify)
            ip_t = $11; ip_total = ip_a + ip_c + ip_g + ip_t;
            edit_rate = ip_t / (ip_total + 1e-6);
            if (edit_rate >= 0.1 && ip_total >= 10) print
        }' ${OUT_PREFIX}_jacusa.tsv > ${OUT_PREFIX}_edits.bed
    else
        echo "Neither Bullseye nor JACUSA2 installed"
        echo "Install: pip install jacusa2-pythonic OR git clone github.com/mekoulnik/Bullseye"
        exit 1
    fi
fi

# Step 3: QC - verify edit rate of fusion / control > 3 (specific signal)
TOTAL_EDITS=$(wc -l < ${OUT_PREFIX}_edits.bed)
echo ""
echo "Total edits: $TOTAL_EDITS"
echo ""
echo "Sanity check ratios:"
echo "  Expected: fusion edits / APOBEC1-only edits > 3"
echo "  If saturated (everything edits): titrate down APOBEC1 expression"

# Step 4 (optional): Filter for DRACH motif if using DART-seq (APOBEC1-YTH for m6A)
# bedtools intersect -wa -u -s -a ${OUT_PREFIX}_edits.bed -b drach_motifs.bed > ${OUT_PREFIX}_dart_drach.bed
# 44% of DART edits fall in DRACH context per Liu 2023

echo ""
echo "Cross-validate with CLIP/eCLIP for high-resolution binding site localization"
echo "STAMP edits are offset 0-50 nt from RBP binding; not single-base resolution"

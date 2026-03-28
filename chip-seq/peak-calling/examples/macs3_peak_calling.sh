#!/bin/bash
# Reference: MACS2 2.2+, MACS3 3.0+, HOMER 4.11+ | Verify API if version differs

OUTPUT_DIR="peaks"
mkdir -p $OUTPUT_DIR

# --- Example 1: MACS3 narrow peaks (TF or H3K4me3) ---
macs3 callpeak \
    -t chip.sorted.bam \
    -c input.sorted.bam \
    -f BAM \
    -g hs \
    -n experiment_narrow \
    --outdir $OUTPUT_DIR \
    -q 0.05 \
    -B --SPMR  # q=0.05 is MACS default FDR; use 0.01 for stricter, 0.1 for exploratory.

# --- Example 2: MACS3 broad peaks (H3K27me3, H3K36me3) ---
macs3 callpeak \
    -t chip.sorted.bam \
    -c input.sorted.bam \
    -f BAM \
    -g hs \
    -n experiment_broad \
    --outdir $OUTPUT_DIR \
    --broad \
    --broad-cutoff 0.1 \
    -B --SPMR  # broad-cutoff=0.1 is MACS default for linking subpeaks; use 0.05 for stricter.

# --- Example 3: tagAlign input, single-chromosome, --nomodel ---
# tagAlign (BED6) is common from ENCODE. Use -f BED.
# For single-chromosome or low-read-count data, skip model building.
# --extsize 150 is typical for H3K4me3; adjust per mark type.
macs3 callpeak \
    -t treatment.tagAlign.gz \
    -c control.tagAlign.gz \
    -f BED \
    -g 46700000 \
    -n chr21_h3k4me3 \
    --outdir $OUTPUT_DIR \
    --nomodel \
    --extsize 150 \
    -q 0.05

# --- Example 4: HOMER peak calling (narrow, TF/H3K4me3) ---
makeTagDirectory chip_tags/ chip.sorted.bam
makeTagDirectory input_tags/ input.sorted.bam
findPeaks chip_tags/ -style factor -i input_tags/ -gsize 2.7e9 -o $OUTPUT_DIR/homer_peaks.txt
pos2bed.pl $OUTPUT_DIR/homer_peaks.txt > $OUTPUT_DIR/homer_peaks.bed

# --- Example 5: HOMER with tagAlign input and custom genome size ---
makeTagDirectory chip_tags_chr21/ treatment.tagAlign.gz -format bed
makeTagDirectory input_tags_chr21/ control.tagAlign.gz -format bed
findPeaks chip_tags_chr21/ -style factor -i input_tags_chr21/ -gsize 46700000 -o $OUTPUT_DIR/homer_chr21.txt
pos2bed.pl $OUTPUT_DIR/homer_chr21.txt > $OUTPUT_DIR/homer_chr21.bed

# --- Example 6: Multi-caller consensus (MACS3 + HOMER intersection) ---
# Intersect MACS3 and HOMER peaks within 500bp for high-confidence set.
bedtools window \
    -a $OUTPUT_DIR/chr21_h3k4me3_peaks.narrowPeak \
    -b $OUTPUT_DIR/homer_chr21.bed \
    -w 500 \
    | cut -f1-5 | sort -k1,1 -k2,2n | uniq > $OUTPUT_DIR/consensus_peaks.bed

# --- Verify results ---
for f in ${OUTPUT_DIR}/*_peaks.narrowPeak ${OUTPUT_DIR}/*_peaks.broadPeak ${OUTPUT_DIR}/homer_*.bed ${OUTPUT_DIR}/consensus_peaks.bed; do
    [ -f "$f" ] && echo "$(basename $f): $(wc -l < $f) peaks"
done

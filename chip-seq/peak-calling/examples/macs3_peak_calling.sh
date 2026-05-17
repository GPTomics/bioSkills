#!/bin/bash
# Reference: MACS2 2.2.9+, MACS3 3.0.4+, HOMER 4.11+, phantompeakqualtools 1.2.2+, IDR 2.0.4+, bedtools 2.31+ | Verify API if version differs
# ChIP-seq peak calling: ENCODE TF pipeline (MACS2 + IDR), ENCODE histone pipeline
# (MACS2 broad + naive overlap), and HOMER alternatives. Assumes deduplicated BAMs
# already filtered with `samtools view -F 1804 -q 30` per ENCODE convention.

set -euo pipefail
OUTPUT_DIR="peaks"
mkdir -p $OUTPUT_DIR

# Effective genome size for hg38, 100bp reads (deepTools tabulated value);
# use 2.913e9 for 50bp, 2.747e9 for 75bp, 2.620e9 for 150bp.
GSIZE=2.701e9


# === ENCODE TF pipeline ===

# 1. Estimate fragment length via cross-correlation (also produces NSC/RSC for QC)
Rscript run_spp.R -c=chip.bam -savp=qc/chip_cc.pdf -out=qc/chip_cc.txt
# Cross-correlation peak; ENCODE pattern reads this from $3 column-1 (estFragLen)
FRAGLEN=$(cut -f3 qc/chip_cc.txt | cut -d',' -f1)

# 2. Per-replicate MACS2 at loose -p 1e-2 threshold (IDR tightens downstream).
# --nomodel --shift 0 --extsize $FRAGLEN is the ENCODE convention regardless of
# whether model-building would succeed; ensures consistency across replicates
# and matches signal-track generation.
for rep in rep1 rep2; do
    macs2 callpeak \
        -t ${rep}.bam -c input.bam \
        -f BAM -g $GSIZE -n $rep \
        --outdir $OUTPUT_DIR/$rep/ \
        --nomodel --shift 0 --extsize $FRAGLEN \
        --keep-dup all \
        -B --SPMR \
        -p 1e-2
done

# 3. Pseudoreplicates: split each rep BAM in half (different seeds!)
samtools view -b -h -s 1.5 rep1.bam > rep1.psr1.bam
samtools view -b -h -s 2.5 rep1.bam > rep1.psr2.bam
samtools view -b -h -s 1.7 rep2.bam > rep2.psr1.bam
samtools view -b -h -s 2.7 rep2.bam > rep2.psr2.bam

# 4. Call peaks on each pseudoreplicate (same parameters)
for psr in rep1.psr1 rep1.psr2 rep2.psr1 rep2.psr2; do
    macs2 callpeak -t ${psr}.bam -c input.bam \
        -f BAM -g $GSIZE -n $psr --outdir $OUTPUT_DIR/$psr/ \
        --nomodel --shift 0 --extsize $FRAGLEN \
        --keep-dup all -p 1e-2
done

# 5. IDR on true replicates (threshold 0.05) -- ENCODE TF standard.
# Sort by column 8 (-log10 p-value); column 7 (signalValue) breaks if MACS
# pile-up scaling differs between libraries.
sort -k8,8nr $OUTPUT_DIR/rep1/rep1_peaks.narrowPeak > $OUTPUT_DIR/rep1.sorted
sort -k8,8nr $OUTPUT_DIR/rep2/rep2_peaks.narrowPeak > $OUTPUT_DIR/rep2.sorted

idr --samples $OUTPUT_DIR/rep1.sorted $OUTPUT_DIR/rep2.sorted \
    --input-file-type narrowPeak --rank p.value \
    --idr-threshold 0.05 \
    --output-file $OUTPUT_DIR/true_reps.idr \
    --plot --log-output-file $OUTPUT_DIR/idr.log

# 6. IDR on pseudoreplicates (looser threshold 0.10) -- per-rep self-consistency
for rep in rep1 rep2; do
    sort -k8,8nr $OUTPUT_DIR/${rep}.psr1/${rep}.psr1_peaks.narrowPeak > $OUTPUT_DIR/${rep}.psr1.sorted
    sort -k8,8nr $OUTPUT_DIR/${rep}.psr2/${rep}.psr2_peaks.narrowPeak > $OUTPUT_DIR/${rep}.psr2.sorted
    idr --samples $OUTPUT_DIR/${rep}.psr1.sorted $OUTPUT_DIR/${rep}.psr2.sorted \
        --input-file-type narrowPeak --rank p.value \
        --idr-threshold 0.10 \
        --output-file $OUTPUT_DIR/${rep}.psr.idr
done

# 7. ENCODE Nself/Nt consistency check (library passes if both ratios <= 2)
Nt=$(awk '$5 >= 540' $OUTPUT_DIR/true_reps.idr | wc -l)  # -125*log2(0.05) ~ 540
N1=$(awk '$5 >= 415' $OUTPUT_DIR/rep1.psr.idr | wc -l)  # -125*log2(0.10) ~ 415
N2=$(awk '$5 >= 415' $OUTPUT_DIR/rep2.psr.idr | wc -l)
echo "Nt=$Nt; N1self=$N1; N2self=$N2"
echo "Self ratio: $(echo "scale=2; if($N2>$N1) $N2/$N1 else $N1/$N2" | bc)"
echo "Rescue ratio: max(Nt,Nself)/min(Nt,Nself); ENCODE rule: both <= 2"


# === ENCODE histone pipeline (broad marks) ===

# H3K27me3, H3K9me3, H3K36me3, H4K20me3: broad mode + naive overlap (NOT IDR;
# IDR's high-vs-low-rank assumption breaks for histone signal dynamic range).
# --broad-cutoff 0.1 controls how MACS stitches sub-peaks within broad regions.
for rep in rep1 rep2 rep3; do
    macs2 callpeak \
        -t ${rep}.bam -c input.bam \
        -f BAM -g $GSIZE -n ${rep}_broad \
        --outdir $OUTPUT_DIR/${rep}_broad/ \
        --broad --broad-cutoff 0.1 \
        --nomodel --shift 0 --extsize $FRAGLEN \
        --keep-dup all -B --SPMR -p 1e-2
done

# Naive overlap: peak appears in >= 2 of 3 replicates with >= 40% reciprocal
# overlap (ENCODE default; commonly misquoted as 50%).
bedtools intersect -a $OUTPUT_DIR/rep1_broad/rep1_broad_peaks.broadPeak \
    -b $OUTPUT_DIR/rep2_broad/rep2_broad_peaks.broadPeak \
    -f 0.40 -r -u > $OUTPUT_DIR/overlap_12.bed
bedtools intersect -a $OUTPUT_DIR/overlap_12.bed \
    -b $OUTPUT_DIR/rep3_broad/rep3_broad_peaks.broadPeak \
    -f 0.40 -r -u > $OUTPUT_DIR/naive_overlap.bed


# === Subset / low-depth fallback (no IDR, single sample) ===

# Single chromosome or pilot data: use tighter -q 0.05 + numeric -g.
# 46.7M is the mappable portion of human chr21 (deepTools approximation).
# --extsize 147 = nucleosome core particle (Buenrostro 2013); biologically
# grounded for nucleosome-proximal marks (H3K4me3, H3K27ac). For TFs use
# the cross-correlation-derived fragment length.
macs3 callpeak \
    -t chr21_chip.tagAlign.gz \
    -c chr21_input.tagAlign.gz \
    -f BED -g 46700000 \
    -n chr21_h3k4me3 \
    --outdir $OUTPUT_DIR/subset/ \
    --nomodel --extsize 147 \
    -q 0.05


# === HOMER alternative (use -style histone for ALL histone marks) ===

# Omnipeak 2025 benchmark: -style histone outperforms -style factor for
# H3K4me3 / H3K27ac because variable-width stitching captures nucleosome-
# adjacent enrichment better than fixed-width factor mode.
makeTagDirectory chip_tags/ chip.bam
makeTagDirectory input_tags/ input.bam

# TF: -style factor (fixed width, three filters: -F 4 -L 4 -C 2)
findPeaks chip_tags/ -style factor -i input_tags/ \
    -gsize $GSIZE -o $OUTPUT_DIR/homer_factor.txt

# Histone mark (any): -style histone (variable width, -L 0)
findPeaks chip_tags/ -style histone -i input_tags/ \
    -gsize $GSIZE -o $OUTPUT_DIR/homer_histone.txt

pos2bed.pl $OUTPUT_DIR/homer_histone.txt > $OUTPUT_DIR/homer_histone.bed


# === Blacklist filter (ENCODE v2; Amemiya 2019) ===

# ENCODE blacklist v2 catches repeat-driven artifacts but NOT hyper-ChIPable
# transcribed genes (Teytelman 2013). For rigorous claims at rRNA/tRNA/histone
# clusters/mitochondrial DNA, build a custom cell-type-specific blacklist from
# the top 1% of input signal.
BLACKLIST_URL="https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz"
wget -nc "$BLACKLIST_URL" || { echo "ERROR: blacklist download failed; check network or URL" >&2; exit 1; }
gunzip -kf hg38-blacklist.v2.bed.gz
bedtools intersect -v -a $OUTPUT_DIR/naive_overlap.bed -b hg38-blacklist.v2.bed \
    > $OUTPUT_DIR/naive_overlap.blacklist_filtered.bed


# === Summary ===
echo ""
echo "=== Peak counts ==="
for f in $OUTPUT_DIR/true_reps.idr $OUTPUT_DIR/naive_overlap*.bed \
         $OUTPUT_DIR/subset/chr21_h3k4me3_peaks.narrowPeak \
         $OUTPUT_DIR/homer_histone.bed; do
    [ -f "$f" ] && echo "$(basename $f): $(wc -l < $f)"
done

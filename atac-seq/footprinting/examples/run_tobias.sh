#!/bin/bash
# Reference: TOBIAS 0.16+, samtools 1.19+, bedtools 2.31+ | Verify API if version differs
# TOBIAS three-step ATAC-seq footprinting with bias correction, scoring, and differential bound/unbound calls.
# Includes CTCF aggregate-footprint sanity check (clean V-shape required for valid downstream calls).

set -euo pipefail

COND1_BAM=${1:-cond1.dedup.nochrM.bam}
COND2_BAM=${2:-cond2.dedup.nochrM.bam}
PEAKS=${3:-consensus_peaks.bed}                  # Must be the same peakset for both conditions
GENOME=${4:-hg38.fa}
BLACKLIST=${5:-hg38-blacklist.v2.bed}
MOTIFS=${6:-JASPAR2024_CORE_vertebrates.pfm}
OUTDIR=${7:-tobias_out}
CORES=${8:-16}

mkdir -p $OUTDIR/{cond1,cond2,bindetect,validation}

# Step 1: Bias correction per condition (matched peakset and blacklist)
for cond in cond1 cond2; do
    BAM=$(eval echo \$${cond^^}_BAM)
    echo "=== ATACorrect: $cond ==="
    TOBIAS ATACorrect \
        --bam $BAM --genome $GENOME \
        --peaks $PEAKS --blacklist $BLACKLIST \
        --outdir $OUTDIR/$cond \
        --cores $CORES
done

# Step 2: Per-base footprint scores
for cond in cond1 cond2; do
    echo "=== ScoreBigwig: $cond ==="
    TOBIAS ScoreBigwig \
        --signal $OUTDIR/$cond/*_corrected.bw \
        --regions $PEAKS \
        --output $OUTDIR/$cond/${cond}_footprints.bw \
        --cores $CORES
done

# Step 3: Differential bound/unbound classification across conditions
echo "=== BINDetect (differential) ==="
TOBIAS BINDetect \
    --motifs $MOTIFS \
    --signals $OUTDIR/cond1/cond1_footprints.bw $OUTDIR/cond2/cond2_footprints.bw \
    --genome $GENOME --peaks $PEAKS \
    --outdir $OUTDIR/bindetect \
    --cond_names cond1 cond2 \
    --cores $CORES

# Validation: aggregate footprint at CTCF (gold-standard QC -- expect clean V-shape)
CTCF_BED=$OUTDIR/bindetect/CTCF_*/beds/CTCF_*_bound.bed   # Whatever motif id JASPAR uses
if [ -f $CTCF_BED ]; then
    TOBIAS PlotAggregate \
        --TFBS $CTCF_BED \
        --signals $OUTDIR/cond1/*_corrected.bw $OUTDIR/cond2/*_corrected.bw \
        --output $OUTDIR/validation/ctcf_aggregate.pdf \
        --share_y both --plot_boundaries
    echo "Validation plot: $OUTDIR/validation/ctcf_aggregate.pdf"
    echo "Expect clean V-shape at CTCF; if shallow, depth or bias correction is the issue."
fi

# Top differential TFs ranked by absolute change
echo "=== Top differential TFs ==="
awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) col[$i]=i; next}
            {print $(col["output_prefix"]), $(col["cond1_cond2_change"]), $(col["cond1_cond2_pvalue"])}' \
    $OUTDIR/bindetect/bindetect_results.txt | \
    sort -k3,3g | head -20

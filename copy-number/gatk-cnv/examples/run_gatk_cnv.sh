#!/usr/bin/env bash
# Reference: GATK 4.5+ | Verify API if version differs
#
# GATK somatic CNV workflow for a tumor-normal pair, including the AnnotateIntervals
# step (explicit GC-bias correction) that is commonly skipped. Output is RELATIVE
# copy-ratio segments with +/-/0 calls and minor-allele fraction -- not integer
# allele-specific copy number, purity, or ploidy. For those, escalate to an
# allele-specific caller (see copy-number/allele-specific-copy-number).

set -euo pipefail

REF=reference.fa                    # needs .fai and .dict alongside
TARGETS=targets.interval_list       # exome targets; for WGS, omit -L and use --bin-length 1000
SNPS=common_biallelic_snps.interval_list
TUMOR_BAM=tumor.bam
NORMAL_BAM=normal.bam
PON=cnv.pon.hdf5                    # built from >= 20 tumor-free, process-matched normals
OUTDIR=gatk_cnv_results
mkdir -p "$OUTDIR/segments" "$OUTDIR/plots"

# 1. Preprocess and annotate intervals. --bin-length 0 uses exome targets as-is.
gatk PreprocessIntervals -R "$REF" -L "$TARGETS" \
    --bin-length 0 --interval-merging-rule OVERLAPPING_ONLY \
    -O "$OUTDIR/preprocessed.interval_list"
gatk AnnotateIntervals -R "$REF" -L "$OUTDIR/preprocessed.interval_list" \
    --interval-merging-rule OVERLAPPING_ONLY -O "$OUTDIR/annotated.tsv"

# 2. Collect read counts for the tumor.
gatk CollectReadCounts -R "$REF" -I "$TUMOR_BAM" \
    -L "$OUTDIR/preprocessed.interval_list" \
    --interval-merging-rule OVERLAPPING_ONLY -O "$OUTDIR/tumor.counts.hdf5"

# 3. Denoise against the panel of normals (tangent normalization). If a known event
# vanishes here but is present in standardizedCR.tsv, the PoN absorbed it -- rebuild
# the PoN larger and tumor-free.
gatk DenoiseReadCounts -I "$OUTDIR/tumor.counts.hdf5" \
    --count-panel-of-normals "$PON" \
    --standardized-copy-ratios "$OUTDIR/tumor.standardizedCR.tsv" \
    --denoised-copy-ratios "$OUTDIR/tumor.denoisedCR.tsv"

# 4. Allelic counts at common SNPs for tumor and matched normal (enables LOH/MAF).
gatk CollectAllelicCounts -R "$REF" -I "$TUMOR_BAM" -L "$SNPS" \
    -O "$OUTDIR/tumor.allelicCounts.tsv"
gatk CollectAllelicCounts -R "$REF" -I "$NORMAL_BAM" -L "$SNPS" \
    -O "$OUTDIR/normal.allelicCounts.tsv"

# 5. Joint segmentation of copy ratio and allele fraction.
gatk ModelSegments --denoised-copy-ratios "$OUTDIR/tumor.denoisedCR.tsv" \
    --allelic-counts "$OUTDIR/tumor.allelicCounts.tsv" \
    --normal-allelic-counts "$OUTDIR/normal.allelicCounts.tsv" \
    --output-prefix tumor -O "$OUTDIR/segments/"

# 6. Call each segment +/-/0 (a simple t-test against the copy-ratio baseline).
gatk CallCopyRatioSegments -I "$OUTDIR/segments/tumor.cr.seg" \
    -O "$OUTDIR/segments/tumor.called.seg"

# 7. Plot denoised ratios and modeled segments.
gatk PlotModeledSegments \
    --denoised-copy-ratios "$OUTDIR/tumor.denoisedCR.tsv" \
    --allelic-counts "$OUTDIR/segments/tumor.hets.tsv" \
    --segments "$OUTDIR/segments/tumor.modelFinal.seg" \
    --sequence-dictionary "${REF%.fa}.dict" \
    --output-prefix tumor -O "$OUTDIR/plots/"

echo "Done. Relative copy-ratio calls: $OUTDIR/segments/tumor.called.seg"

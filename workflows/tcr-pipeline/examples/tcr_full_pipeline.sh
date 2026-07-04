#!/bin/bash
# Reference: MiXCR 4.x VDJtools 1.2.x | Verify API if version differs
# End-to-end immune-repertoire pipeline, fork-aware (bulk vs single-cell, TCR vs BCR).
# MiXCR 4.x is preset-driven and license-gated; the 3.x `mixcr align -s hsa -p rna-seq`
# chain is removed. Bulk diversity is compared ONLY after downsampling to a common depth.
set -euo pipefail

R1=$1
R2=$2
SAMPLE=$3
OUTDIR=${4:-'tcr_results'}
PRESET=${5:-'generic-amplicon'}   # match to the exact kit; the wrong preset silently corrupts CDR3
SPECIES=${6:-'hsa'}               # required for generic presets: hsa (human), mmu (mouse)
RECEPTOR=${7:-'TCR'}              # TCR -> exact clonotypes + VDJtools; BCR -> Immcantation
CHAIN=${8:-'TRB'}                 # TRB/TRA for TCR, IGH for BCR
DOWNSAMPLE_TO=${9:-50000}         # common depth for cross-sample diversity; set to smallest sample

mkdir -p "${OUTDIR}"/{mixcr,vdjtools,airr,plots}

# Stage 0: license (academic is free). Skips if already activated or MI_LICENSE_FILE is set.
mixcr activate-license 2>/dev/null || echo 'License already active or MI_LICENSE_FILE set'

echo '=== Stage 1: MiXCR analyze (align -> refine -> assemble in one command) ==='
# From 4.7, presets without an intrinsic assembling feature need --assemble-clonotypes-by CDR3.
mixcr analyze "${PRESET}" \
    --species "${SPECIES}" \
    -f \
    "${R1}" "${R2}" \
    "${OUTDIR}/mixcr/${SAMPLE}"

echo '=== Stage 1 QC: alignment rate + chain usage (catch wrong preset / contamination) ==='
mixcr qc "${OUTDIR}/mixcr/${SAMPLE}.clns"
mixcr exportQc align "${OUTDIR}/mixcr/${SAMPLE}.clns" "${OUTDIR}/plots/qc_align.pdf"
mixcr exportQc chainUsage "${OUTDIR}/mixcr/${SAMPLE}.clns" "${OUTDIR}/plots/qc_chains.pdf"

echo '=== Stage 2: export (VDJtools table for bulk TCR, AIRR TSV for BCR/single-cell) ==='
mixcr exportClones -c "${CHAIN}" \
    "${OUTDIR}/mixcr/${SAMPLE}.clns" \
    "${OUTDIR}/mixcr/${SAMPLE}.clones_${CHAIN}.tsv"
mixcr exportAirr \
    "${OUTDIR}/mixcr/${SAMPLE}.clns" \
    "${OUTDIR}/airr/${SAMPLE}.airr.tsv"

if [ "${RECEPTOR}" = 'BCR' ]; then
    echo '=== BCR fork: exact clonotypes are WRONG (somatic hypermutation) ==='
    echo 'Hand the AIRR TSV to Immcantation (tcr-bcr-analysis/immcantation-analysis):'
    echo '  distToNearest -> findThreshold -> hierarchicalClones'
    echo '  -> CreateGermlines.py --cloned -> observedMutations -> dowser getTrees'
    echo "AIRR: ${OUTDIR}/airr/${SAMPLE}.airr.tsv"
    exit 0
fi

echo '=== Stage 3t (bulk TCR): convert, then DOWNSAMPLE before any diversity/overlap ==='
# Diversity, clonality and Jaccard overlap all grow with depth; comparing raw values across
# unequal-depth samples measures depth, not biology. Equalize depth first.
vdjtools Convert -S mixcr \
    "${OUTDIR}/mixcr/${SAMPLE}.clones_${CHAIN}.tsv" \
    "${OUTDIR}/vdjtools/"

vdjtools DownSample -x "${DOWNSAMPLE_TO}" \
    "${OUTDIR}/vdjtools/${SAMPLE}.clones_${CHAIN}.tsv" \
    "${OUTDIR}/vdjtools/ds_"

# Diversity on the depth-normalized sample; report the resampled table for cross-sample claims.
vdjtools CalcDiversityStats \
    "${OUTDIR}/vdjtools/ds_${SAMPLE}.clones_${CHAIN}.tsv" \
    "${OUTDIR}/vdjtools/diversity"

echo '=== Stage 4: visualization ==='
vdjtools PlotFancyVJUsage \
    "${OUTDIR}/vdjtools/ds_${SAMPLE}.clones_${CHAIN}.tsv" \
    "${OUTDIR}/plots/vj_usage"
vdjtools PlotFancySpectratype \
    "${OUTDIR}/vdjtools/ds_${SAMPLE}.clones_${CHAIN}.tsv" \
    "${OUTDIR}/plots/spectratype"

echo '=== Pipeline complete ==='
echo "Clonotypes: ${OUTDIR}/mixcr/${SAMPLE}.clones_${CHAIN}.tsv"
echo "AIRR:       ${OUTDIR}/airr/${SAMPLE}.airr.tsv"
echo "Diversity:  ${OUTDIR}/vdjtools/diversity.diversity.strict.txt"
echo "Plots:      ${OUTDIR}/plots/"
echo 'For multi-sample cohorts: run DownSample with a metadata file (-m) so all samples share one depth.'

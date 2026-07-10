#!/bin/bash
# Reference: MiXCR 4.x | Verify API if version differs
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
RECEPTOR=${7:-'TCR'}              # TCR -> exact clonotypes + MiXCR postanalysis; BCR -> Immcantation
CHAIN=${8:-'TRB'}                 # TRB/TRA for TCR, IGH for BCR
DOWNSAMPLE_TO=${9:-50000}         # common depth for cross-sample diversity; set near the cohort lower quartile and EXCLUDE (do not normalize down to) samples far below it
# The weight/downsample unit MUST match the library chemistry. Weighting a UMI library by reads
# re-introduces the exact PCR bias the UMIs removed, inflating clonality and deflating diversity.
# Auto-detect from the preset; override with $11 for kits whose names do not encode the chemistry.
case "${PRESET}" in
    *umi*)                 COUNT_UNIT='umi' ;;
    *sc*|*10x*|*rhapsody*) COUNT_UNIT='cell' ;;
    *)                     COUNT_UNIT='read' ;;
esac
COUNT_UNIT=${11:-${COUNT_UNIT}}
# DOWNSAMPLE_TO is in COUNT_UNIT units: molecule/cell counts run ~1-2 orders below read counts,
# so a read-scale default will drop every sample on a UMI or single-cell run. Retune per chemistry.
# generic-amplicon REQUIRES material type + both boundary mixins (errors without them). Set to ''
# for self-contained presets (Takara/BD/10x). 5'RACE libraries use --rigid-left-alignment-boundary.
MIXINS=${10:-'--rna --floating-left-alignment-boundary --floating-right-alignment-boundary C'}

mkdir -p "${OUTDIR}"/{mixcr,postanalysis,airr,plots}

# Stage 0: license (academic is free). Skips if already activated or MI_LICENSE_FILE is set.
mixcr activate-license 2>/dev/null || echo 'License already active or MI_LICENSE_FILE set'

echo '=== Stage 1: MiXCR analyze (align -> refine -> assemble in one command) ==='
# From 4.7, presets without an intrinsic assembling feature need --assemble-clonotypes-by CDR3.
mixcr analyze "${PRESET}" \
    --species "${SPECIES}" \
    ${MIXINS} \
    -f \
    "${R1}" "${R2}" \
    "${OUTDIR}/mixcr/${SAMPLE}"

echo '=== Stage 1 QC: alignment rate + chain usage (catch wrong preset / contamination) ==='
mixcr qc "${OUTDIR}/mixcr/${SAMPLE}.clns"
mixcr exportQc align "${OUTDIR}/mixcr/${SAMPLE}.clns" "${OUTDIR}/plots/qc_align.pdf"
mixcr exportQc chainUsage "${OUTDIR}/mixcr/${SAMPLE}.clns" "${OUTDIR}/plots/qc_chains.pdf"

echo '=== Stage 2: export (MiXCR clones table for bulk TCR, AIRR TSV for BCR/single-cell) ==='
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

echo '=== Stage 3t (bulk TCR): DOWNSAMPLE before any diversity/overlap ==='
# Diversity, clonality and Jaccard overlap all grow with depth; comparing raw values across
# unequal-depth samples measures depth, not biology. Equalize depth first.
# VDJtools is unmaintained for MiXCR 4.x: its parser fails on raw 4.x exportClones output
# ("Unable to parse clonotype string"). Use the native postanalysis, which consumes .clns directly
# and applies the SAME downsample-first semantics (--default-downsampling runs before every metric).
mixcr postanalysis individual \
    --default-downsampling "count-${COUNT_UNIT}-fixed-${DOWNSAMPLE_TO}" \
    --default-weight-function "${COUNT_UNIT}" \
    --only-productive \
    --chains "${CHAIN}" \
    --tables "${OUTDIR}/postanalysis/tables.tsv" \
    "${OUTDIR}/mixcr/${SAMPLE}.clns" \
    "${OUTDIR}/postanalysis/pa.json.gz"

echo '=== Stage 4: visualization ==='
mixcr exportPlots diversity \
    "${OUTDIR}/postanalysis/pa.json.gz" "${OUTDIR}/plots/diversity.pdf"
mixcr exportPlots vjUsage \
    "${OUTDIR}/postanalysis/pa.json.gz" "${OUTDIR}/plots/vj_usage.pdf"

echo '=== Pipeline complete ==='
echo "Clonotypes: ${OUTDIR}/mixcr/${SAMPLE}.clones_${CHAIN}.tsv"
echo "AIRR:       ${OUTDIR}/airr/${SAMPLE}.airr.tsv"
echo "Diversity:  ${OUTDIR}/postanalysis/tables.tsv (+ pa.json.gz)"
echo "Plots:      ${OUTDIR}/plots/"
echo 'For multi-sample cohorts: pass all .clns files (or their directory) to one postanalysis run with --metadata/--group so every sample shares one downsampling depth.'

#!/bin/bash
# Reference: MiXCR 4.x | Verify API if version differs
# MiXCR 4.x TCR/BCR clonotype workflow: preset-driven analyze -> QC -> export.
# The preset encodes the chemistry (material, boundaries, UMI/cell barcodes) and IS
# ~90% of correctness; a mismatched preset fails SILENTLY. Audit it before trusting output.
set -euo pipefail

R1=${1:?usage: mixcr_workflow.sh R1.fastq.gz R2.fastq.gz [preset] [out_prefix] [species] [chain]}
R2=${2:?need R2}
PRESET=${3:-generic-amplicon}          # e.g. 10x-sc-xcr-vdj, takara-human-rna-tcr-umi-smarter-v2, rna-seq
OUT=${4:-sample}
SPECIES=${5:-hsa}                       # hsa|mmu|... required for generic-* presets
CHAIN=${6:-TRB}                         # TRA|TRB|TRG|TRD|IGH|IGK|IGL

# MiXCR 4.x refuses to run unlicensed. Academic use is free (platforma.bio/getlicense).
# Set MI_LICENSE_FILE (or MI_LICENSE), or run: mixcr activate-license
if ! mixcr --version >/dev/null 2>&1; then
    echo 'mixcr not on PATH' >&2; exit 1
fi

# Audit exactly what the preset does BEFORE running (resolved parameter YAML).
mixcr exportPreset --preset-name "$PRESET" "${OUT}.preset.yaml"

# One-command pipeline. --assemble-clonotypes-by is required (4.7+) when the preset
# has no intrinsic assembling feature; CDR3 is the robust default on short reads.
# --species matters for generic-* presets. -f overwrites prior outputs.
mixcr analyze "$PRESET" \
    --species "$SPECIES" \
    --assemble-clonotypes-by CDR3 \
    -f \
    "$R1" "$R2" \
    "$OUT"

# QC: alignment rate + chain composition. Low align rate or off-target chains means
# wrong preset/species/boundaries or cross-contamination -- caught only here.
mixcr qc "${OUT}.clns"
mixcr exportQc align "${OUT}.clns" "${OUT}.qc_align.pdf"
mixcr exportQc chainUsage "${OUT}.clns" "${OUT}.qc_chainUsage.pdf"

# Native MiXCR export. readCount for non-UMI bulk; uniqueMoleculeCount is the correct
# abundance on UMI libraries (reporting reads there re-adds PCR bias). D-gene is
# near-unassignable in TRB/IGH, so it is deliberately NOT exported as a key field.
mixcr exportClones -c "$CHAIN" \
    -cloneId -readCount -readFraction -uniqueMoleculeCount \
    -nSeqCDR3 -aaSeqCDR3 -bestVGene -bestJGene -allVHitsWithScore -isProductive VRegion \
    "${OUT}.clns" "${OUT}.clones_${CHAIN}.tsv"

# AIRR rearrangement TSV for Immcantation / scirpy / any AIRR tool. Use the dedicated
# exportAirr command -- downstream tools expect AIRR headers, not renamed native ones.
mixcr exportAirr "${OUT}.clns" "${OUT}.airr.tsv"

echo "Done. Preset audited in ${OUT}.preset.yaml; clonotypes in ${OUT}.clones_${CHAIN}.tsv and ${OUT}.airr.tsv"

#!/bin/bash
# Reference: isoseq 4.3+, pbmm2 1.13+, pigeon 1.2+ | Verify API if version differs
# PacBio Iso-Seq / Kinnex full-length isoform pipeline. Classification + filter IS the
# analysis: a novel isoform is an artifact until orthogonal support clears it.
set -euo pipefail

HIFI_BAM=${1:?Usage: $0 <hifi_reads.bam> <primers.fa> <ref.fa> <annotation.gtf> [is_kinnex] [cage.bed] [polyA.txt]}
PRIMERS=${2:?cDNA primers FASTA required}
REFERENCE=${3:?reference genome required}
ANNOTATION=${4:?reference annotation GTF required}
IS_KINNEX=${5:-yes}        # Kinnex/MAS-seq needs skera deconcatenation first
CAGE=${6:-}                # CAGE refTSS BED for 5' TSS validation (catches 5' degradation)
POLYA=${7:-}               # poly-A motif list for 3' validation (catches intra-priming)

# 0. Kinnex only: deconcatenate the MAS array into segmented reads before lima.
if [ "$IS_KINNEX" = yes ]; then
    skera split "$HIFI_BAM" mas_adapters.fasta segmented.bam
    INPUT=segmented.bam
else
    INPUT="$HIFI_BAM"
fi

# 1. primers -> 2. FLNC (full-length non-chimeric, poly-A required) -> 3. map -> 4. collapse
lima "$INPUT" "$PRIMERS" fl.bam --isoseq --peek-guess
isoseq refine fl.5p--3p.bam "$PRIMERS" flnc.bam --require-polya
pbmm2 align --preset ISOSEQ --sort "$REFERENCE" flnc.bam mapped.bam
isoseq collapse --do-not-collapse-extra-5exons mapped.bam flnc.bam collapsed.gff
#   collapsed.flnc_count.txt = FLNC molecules per isoform = the real depth metric

# 5. classify + filter with pigeon (consumes the collapsed.sorted.gff, NOT a BAM)
pigeon prepare collapsed.gff                       # sorts the transcript GFF
pigeon prepare "$ANNOTATION" "$REFERENCE"          # sorts the annotation -> annotation.sorted.gtf, indexes genome
CAGE_ARG=""; [ -n "$CAGE" ] && CAGE_ARG="--cage-peak $CAGE"
POLYA_ARG=""; [ -n "$POLYA" ] && POLYA_ARG="--poly-a $POLYA"
pigeon classify collapsed.sorted.gff "${ANNOTATION%.gtf}.sorted.gtf" "$REFERENCE" \
    --fl collapsed.flnc_count.txt $CAGE_ARG $POLYA_ARG
pigeon filter collapsed_classification.txt --isoforms collapsed.sorted.gff

# Saturation excluding singletons (1-FLNC isoforms are the unreliable bucket)
pigeon report --exclude-singletons \
    collapsed_classification.filtered_lite_classification.txt saturation.txt
echo 'Done. Inspect the FSM/ISM ratio (high ISM = 5-prime degradation, not novel biology).'

#!/usr/bin/env bash
# Reference: MACS3 3.0+, samtools 1.19+, bedtools 2.31+ | Verify with `macs3 --version`, `samtools --version`, `bedtools --version` if installed releases differ.
# Broad m6A peak calling with MACS3 — second-opinion or viral / kb-scale peaks.
# CRITICAL: --keep-dup all is mandatory for MeRIP (default --keep-dup 1 destroys signal).

set -euo pipefail

IP_BAMS=(aligned/IP_rep1_Aligned.sortedByCoord.out.bam
         aligned/IP_rep2_Aligned.sortedByCoord.out.bam
         aligned/IP_rep3_Aligned.sortedByCoord.out.bam)

INPUT_BAMS=(aligned/Input_rep1_Aligned.sortedByCoord.out.bam
            aligned/Input_rep2_Aligned.sortedByCoord.out.bam
            aligned/Input_rep3_Aligned.sortedByCoord.out.bam)

OUTDIR='macs3_output'
NAME='m6a_run1'
GSIZE='hs'                # hs for human (~2.7e9); mm for mouse (~1.87e9)
EXTSIZE=150               # MeRIP fragment length convention
BROAD_CUTOFF=0.1          # MACS2/3 default for broad mode
QVALUE=0.05

mkdir -p "${OUTDIR}"

macs3 callpeak \
    --treatment "${IP_BAMS[@]}" \
    --control "${INPUT_BAMS[@]}" \
    --format BAMPE \
    --gsize "${GSIZE}" \
    --nomodel \
    --extsize "${EXTSIZE}" \
    --keep-dup all \
    --broad \
    --broad-cutoff "${BROAD_CUTOFF}" \
    --qvalue "${QVALUE}" \
    --outdir "${OUTDIR}" \
    --name "${NAME}"

narrowpeak="${OUTDIR}/${NAME}_peaks.narrowPeak"
broadpeak="${OUTDIR}/${NAME}_peaks.broadPeak"

if [ -f "${narrowpeak}" ]; then
    echo "MACS3 narrow peaks: $(wc -l < ${narrowpeak})"
fi

if [ -f "${broadpeak}" ]; then
    echo "MACS3 broad peaks: $(wc -l < ${broadpeak})"
fi

# Intersect with exomePeak2 peaks (if available) for 2-tool consensus.
EXOMEPEAK_BED='exomepeak2_output/m6a_run1/peaks.bed'

if [ -f "${EXOMEPEAK_BED}" ] && [ -f "${broadpeak}" ]; then
    bedtools intersect \
        -a "${EXOMEPEAK_BED}" \
        -b "${broadpeak}" \
        -wa -u -f 0.5 \
        > "${OUTDIR}/${NAME}_consensus_exomepeak_macs3.bed"

    echo "exomePeak2 + MACS3 broad consensus: $(wc -l < ${OUTDIR}/${NAME}_consensus_exomepeak_macs3.bed)"
fi

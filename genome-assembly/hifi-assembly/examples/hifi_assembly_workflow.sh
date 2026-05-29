#!/bin/bash
# Reference: hifiasm 0.25.0+, yak 0.1+, gfatools 0.5+ | Verify API if version differs
# Phased diploid HiFi assembly with hifiasm: pick the phasing mode by available data,
# extract haplotypes from GFA, and flag that phasing QC needs k-mer/trio metrics (see assembly-qc).

set -euo pipefail

SAMPLE="sample"
HIFI_READS="reads.hifi.fastq.gz"
THREADS=32

# Phasing inputs - leave empty for HiFi-only (partial phasing).
HIC_R1=""                 # Hi-C R1; if set with HIC_R2 -> global Hi-C phasing (no parents needed)
HIC_R2=""                 # Hi-C R2
PAT_READS=""              # paternal Illumina (trio); if set with MAT_READS -> gold-standard phasing
MAT_READS=""              # maternal Illumina (trio)

# Sample-type knobs - the load-bearing decisions, NOT cosmetic defaults.
PURGE_LEVEL=3             # hifiasm HiFi-only/Hi-C default; set 0 for inbred/doubled-haploid/mole or real
                          # segmental duplications get DELETED as if they were duplicate haplotigs.
                          # NOTE: trio mode intentionally ignores this (hifiasm defaults trio to -l0;
                          # parental binning, not purging, separates the haplotypes).
YAK_K=31                  # yak/hifiasm trio k-mer size; do NOT use 30 (that is verkko/Merqury hap-mers)
HOM_COV=0                 # 0 = let hifiasm estimate; set to the k-mer-histogram peak if hap sizes come
                          # out unbalanced (an unbalanced ratio is a coverage-estimate alarm, not biology)

gfa_to_fasta () {
    # Extract contig sequences from GFA S (segment) lines. hifiasm emits GFA, not FASTA.
    gfatools gfa2fa "$1" > "$2" 2>/dev/null || awk '/^S/{print ">"$2"\n"$3}' "$1" > "$2"
}

echo '=== Step 1: pick phasing mode from available data ==='
# -l (purge) is added per-branch: trio omits it (hifiasm defaults trio to -l0).
HIFIASM_ARGS=(-o "${SAMPLE}" -t "${THREADS}")
[ "${HOM_COV}" -gt 0 ] && HIFIASM_ARGS+=(--hom-cov "${HOM_COV}")

if [ -n "${PAT_READS}" ] && [ -n "${MAT_READS}" ]; then
    echo 'trio mode (gold standard) -> outputs use the .dip. prefix'
    yak count -k"${YAK_K}" -b37 -t16 -o pat.yak ${PAT_READS}
    yak count -k"${YAK_K}" -b37 -t16 -o mat.yak ${MAT_READS}
    HIFIASM_ARGS+=(-1 pat.yak -2 mat.yak)   # no -l: trio defaults to -l0 (binning, not purging)
    PREFIX="${SAMPLE}.dip"
elif [ -n "${HIC_R1}" ] && [ -n "${HIC_R2}" ]; then
    echo 'Hi-C mode (global phasing, no parents) -> outputs use the .bp. prefix'
    HIFIASM_ARGS+=(--h1 "${HIC_R1}" --h2 "${HIC_R2}" -l "${PURGE_LEVEL}")
    PREFIX="${SAMPLE}.bp"
else
    echo 'HiFi-only mode -> PARTIAL phasing (switch errors between blocks); .bp. prefix'
    HIFIASM_ARGS+=(-l "${PURGE_LEVEL}")
    PREFIX="${SAMPLE}.bp"
fi

echo '=== Step 2: assemble ==='
hifiasm "${HIFIASM_ARGS[@]}" "${HIFI_READS}"

echo '=== Step 3: extract haplotypes + primary from GFA ==='
gfa_to_fasta "${PREFIX}.hap1.p_ctg.gfa" "${SAMPLE}.hap1.fasta"
gfa_to_fasta "${PREFIX}.hap2.p_ctg.gfa" "${SAMPLE}.hap2.fasta"
# The primary exists only in .bp. (HiFi-only/Hi-C) mode; it is a MOSAIC, not a haplotype.
[ -f "${PREFIX}.p_ctg.gfa" ] && gfa_to_fasta "${PREFIX}.p_ctg.gfa" "${SAMPLE}.primary.fasta"

echo '=== Step 4: next steps (NOT a phasing check) ==='
echo 'Contiguity/completeness (N50, BUSCO) and -- critically -- PHASING QC live in assembly-qc.'
echo 'Switch errors are invisible to N50/BUSCO/QV: validate with Merqury hap-mer blob plots'
echo 'and switch/hamming error against trio or Hi-C hap-mers. No hap-mers = phasing unvalidated.'

echo "Done. hap1=${SAMPLE}.hap1.fasta hap2=${SAMPLE}.hap2.fasta (use phased haps, NOT primary, for allele-aware work)"

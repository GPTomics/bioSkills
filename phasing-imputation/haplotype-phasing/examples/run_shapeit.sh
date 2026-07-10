#!/bin/bash
# Reference: SHAPEIT5 5.1.1, bcftools 1.19+ | Verify API if version differs
#
# Statistical phasing with the SHAPEIT5 common-scaffold-then-rare design:
# phase_common (common variants -> scaffold) -> ligate (stitch overlapping chunks)
# -> phase_rare (place rare variants on the fixed scaffold) -> switch (benchmark vs
# a trio). NOT a spot-run target: it needs the SHAPEIT5 binaries, a build-matched
# genetic map, and minutes-to-hours of compute on real data. Parse-check and verify
# every flag against the installed `phase_common --help` before running.
set -euo pipefail

TARGET=${1:?'usage: run_shapeit.sh target.bcf chr20.b38.gmap.gz chr20'}
MAP=${2:?'genetic map matching the data genome build'}
REGION=${3:?'chromosome or region, e.g. chr20'}
OUT=${OUT:-phasing_out}
THREADS=${THREADS:-16}

COMMON_MAF=0.001   # common/rare split: variants above this build the scaffold; rarer ones are phased onto it. SHAPEIT5 default boundary for the scaffold
# shellcheck disable=SC2034  # surfaced as a tunable decision threshold for the agent, not consumed by a command here
RARE_N_MIN=2000    # below this sample size the phase_rare step adds little (too few rare-allele carriers); phase_common alone suffices

mkdir -p "${OUT}"

# 1. Phase COMMON variants into a scaffold (run per overlapping chunk for big chromosomes)
phase_common \
    --input "${TARGET}" \
    --filter-maf "${COMMON_MAF}" \
    --region "${REGION}" \
    --map "${MAP}" \
    --output "${OUT}/scaffold.${REGION}.bcf" \
    --thread "${THREADS}"

# 2. Ligate per-chunk scaffolds into one chromosome (chunks MUST overlap; skip if single-chunk)
# ligate --input "${OUT}/scaffold_chunks.txt" --output "${OUT}/scaffold.bcf" --thread "${THREADS}" --index

# 3. Phase RARE variants onto the fixed scaffold using the FULL genotypes
phase_rare \
    --input "${TARGET}" \
    --scaffold "${OUT}/scaffold.${REGION}.bcf" \
    --map "${MAP}" \
    --input-region "${REGION}" \
    --scaffold-region "${REGION}" \
    --output "${OUT}/phased.${REGION}.bcf" \
    --thread "${THREADS}"

# 4. Benchmark against a trio truth set (Mendelian phase); reports SER/GER stratified by MAC
# switch --validation trio_truth.bcf --estimation "${OUT}/phased.${REGION}.bcf" \
#        --pedigree family.fam --region "${REGION}" --output "${OUT}/eval" --thread "${THREADS}"

echo "== phased haplotypes in ${OUT}/phased.${REGION}.bcf; report switch error stratified by MAC, not genome-wide =="

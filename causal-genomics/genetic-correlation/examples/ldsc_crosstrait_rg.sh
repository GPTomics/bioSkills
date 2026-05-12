#!/usr/bin/env bash
# Reference: LDSC v1.0.1+, HDL, LAVA | Verify API if version differs
# Cross-trait LDSC pipeline for global genetic correlation.
# Inputs are two raw GWAS sumstats files; outputs rg + intercept + ratio.
#
# Operational rule: cross-trait LDSC intercept absorbs sample overlap WITHOUT
# biasing the rg estimate (Bulik-Sullivan 2015 Nat Genet 47:1236). Report the
# intercept as an overlap diagnostic but do not interpret it as a bias on rg.

set -euo pipefail

TRAIT1_RAW=${1:?trait1 raw sumstats}
TRAIT2_RAW=${2:?trait2 raw sumstats}
TRAIT1_N=${3:?trait1 effective N}
TRAIT2_N=${4:?trait2 effective N}
LDSC_REF=${5:-eur_w_ld_chr/}   # ancestry-matched; EUR default
HM3_SNPS=${6:-w_hm3.snplist}   # HapMap3 SNPs for munge restriction
OUT_PREFIX=${7:-rg_t1_t2}

# Underpower floor: LDSC mean chi-square > 1.02 required for stable rg
# (LDSC tutorial); below this, no method recovers precision.
MIN_MEAN_CHISQ=1.02

# Local-rg trigger: |rg| > 0.5 from global LDSC warrants LAVA local rg per
# locus and (if MR is downstream) CHP-aware sensitivity (CAUSE / LHC-MR).
RG_LOCAL_TRIGGER=0.5

echo '== Step 1: munge trait1 =='
munge_sumstats.py \
    --sumstats "${TRAIT1_RAW}" \
    --N "${TRAIT1_N}" \
    --merge-alleles "${HM3_SNPS}" \
    --out "${OUT_PREFIX}.t1"

echo '== Step 2: munge trait2 =='
munge_sumstats.py \
    --sumstats "${TRAIT2_RAW}" \
    --N "${TRAIT2_N}" \
    --merge-alleles "${HM3_SNPS}" \
    --out "${OUT_PREFIX}.t2"

echo '== Step 3: pre-check mean chi-square (each trait) =='
for tag in t1 t2; do
    munged="${OUT_PREFIX}.${tag}.sumstats.gz"
    mean_chisq=$(zcat "${munged}" | awk 'NR>1 {sum += $5*$5; n++} END {if (n>0) print sum/n; else print 0}')
    echo "trait ${tag} mean chi-square: ${mean_chisq}"
    pass=$(awk -v m="${mean_chisq}" -v thr="${MIN_MEAN_CHISQ}" 'BEGIN {print (m > thr) ? 1 : 0}')
    if [[ "${pass}" == "0" ]]; then
        echo "WARNING: trait ${tag} mean chi-square below ${MIN_MEAN_CHISQ}; LDSC rg will be underpowered."
    fi
done

echo '== Step 4: cross-trait LDSC rg =='
ldsc.py \
    --rg "${OUT_PREFIX}.t1.sumstats.gz,${OUT_PREFIX}.t2.sumstats.gz" \
    --ref-ld-chr "${LDSC_REF}" \
    --w-ld-chr "${LDSC_REF}" \
    --out "${OUT_PREFIX}"

echo '== Step 5: extract rg block =='
grep -A 11 'Summary of Genetic Correlation Results' "${OUT_PREFIX}.log" || true

echo '== Step 6: interpret intercept (overlap diagnostic) =='
gcov_int=$(grep 'gcov_int:' "${OUT_PREFIX}.log" | head -1 | awk '{print $2}')
rg_pt=$(grep -A 11 'Summary of Genetic Correlation Results' "${OUT_PREFIX}.log" | tail -1 | awk '{print $3}')
echo "rg = ${rg_pt}"
echo "gcov_int (sample-overlap proxy; NOT a bias on rg) = ${gcov_int}"

echo '== Step 7: flag LAVA / CHP-aware MR triggers =='
if [[ -n "${rg_pt}" ]]; then
    abs_rg=$(awk -v r="${rg_pt}" 'BEGIN {print (r < 0) ? -r : r}')
    above=$(awk -v a="${abs_rg}" -v t="${RG_LOCAL_TRIGGER}" 'BEGIN {print (a > t) ? 1 : 0}')
    if [[ "${above}" == "1" ]]; then
        echo "ACTION: |rg| > ${RG_LOCAL_TRIGGER} -> run LAVA local rg (lava_local_rg.R) AND add CHP-aware MR sensitivity (causal-genomics/pleiotropy-detection)."
    else
        echo "Global rg modest; if biology suggests shared etiology, still run LAVA (cancellation can hide local rg)."
    fi
fi

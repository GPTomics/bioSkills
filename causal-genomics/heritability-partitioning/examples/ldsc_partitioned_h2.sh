#!/bin/bash
# Reference: LDSC v1.0.1+, baselineLD_v2.2, Multi_tissue_chromatin_1000Gv3 | Verify API if version differs
#
# End-to-end LDSC pipeline:
#   1. Munge GWAS summary statistics into LDSC format
#   2. Total h2 (with intercept / ratio / mean chi-square)
#   3. Partitioned h2 across baseline-LD v2.2 functional annotations (Gazal 2017)
#   4. Cell-type prioritization via Finucane 2018 chromatin .ldcts manifest
#   5. Cross-trait genetic correlation (rg) example
#
# Usage:
#   bash ldsc_partitioned_h2.sh <gwas.tsv> <trait_prefix> [samp_prev] [pop_prev]
#
# Prerequisite environment:
#   Python 3 LDSC fork installed; prefer abdenlab/ldsc-python3 (v2.0.0) because
#   belowlab/ldsc v3.0.1 README states the --h2/--rg/--h2-cts CLI is broken (use
#   Docker jtb114/ldsc:latest for the belowlab fallback). See SKILL.md.
#   Reference resources downloaded from https://alkesgroup.broadinstitute.org/LDSCORE/
#     - eur_w_ld_chr/                                          (univariate h2)
#     - 1000G_Phase3_baselineLD_v2.2_ldscores/                 (partitioned h2)
#     - 1000G_Phase3_frq/1000G.EUR.QC.{1..22}.frq              (allele frequencies)
#     - 1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.     (regression weights)
#     - Multi_tissue_chromatin_1000Gv3_ldscores/               (Finucane 2018 cell-type)
#     - w_hm3.snplist                                          (HapMap3 SNP filter for munge)

set -euo pipefail

GWAS_FILE=${1:?Usage: ldsc_partitioned_h2.sh <gwas.tsv> <trait_prefix> [samp_prev] [pop_prev]}
TRAIT=${2:?provide trait_prefix}
SAMP_PREV=${3:-}     # case fraction in GWAS sample; empty for quantitative trait
POP_PREV=${4:-}      # population lifetime prevalence; empty for quantitative trait

REF_DIR=${LDSC_REF:-./ldsc_refs}
EUR_LD=${REF_DIR}/eur_w_ld_chr
BASELINE=${REF_DIR}/1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD.
FRQ=${REF_DIR}/1000G_Phase3_frq/1000G.EUR.QC.
WEIGHTS=${REF_DIR}/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.
CTS_MANIFEST=${REF_DIR}/Multi_tissue_chromatin_1000Gv3_ldscores/Multi_tissue_chromatin.ldcts
HM3_SNPS=${REF_DIR}/w_hm3.snplist

OUT_DIR=${TRAIT}_ldsc
mkdir -p ${OUT_DIR}

# ============================================================
# Step 1: Munge sumstats to LDSC format
# ============================================================
# munge_sumstats.py filters to HapMap3 SNPs, harmonises alleles, computes Z scores
# Required columns: SNP A1 A2 BETA (or Z) SE P N (per-row N or supply --N flag)

munge_sumstats.py \
    --sumstats ${GWAS_FILE} \
    --merge-alleles ${HM3_SNPS} \
    --snp SNP --a1 A1 --a2 A2 \
    --signed-sumstats BETA,0 \
    --p P --N-col N \
    --out ${OUT_DIR}/${TRAIT}

MUNGED=${OUT_DIR}/${TRAIT}.sumstats.gz

# ============================================================
# Step 2: Total h2 (univariate)
# ============================================================
# For case-control, supply --samp-prev and --pop-prev to get liability-scale h2

H2_FLAGS="--h2 ${MUNGED} --ref-ld-chr ${EUR_LD}/ --w-ld-chr ${EUR_LD}/"

if [[ -n "${SAMP_PREV}" && -n "${POP_PREV}" ]]; then
    H2_FLAGS="${H2_FLAGS} --samp-prev ${SAMP_PREV} --pop-prev ${POP_PREV}"
    echo "Computing liability-scale h2 (samp_prev=${SAMP_PREV} pop_prev=${POP_PREV})"
fi

ldsc.py ${H2_FLAGS} --out ${OUT_DIR}/${TRAIT}_h2

# Inspect the .log:
#   Total Observed scale h2 (or Liability scale)
#   Lambda GC
#   Mean Chi^2
#   Intercept
#   Ratio = (intercept - 1) / (mean_chi2 - 1)
# Ratio < 0.2 -> mostly polygenic, h2 trustworthy

# ============================================================
# Step 3: Partitioned h2 with baseline-LD v2.2 (Gazal 2017)
# ============================================================

ldsc.py \
    --h2 ${MUNGED} \
    --ref-ld-chr ${BASELINE} \
    --frqfile-chr ${FRQ} \
    --w-ld-chr ${WEIGHTS} \
    --overlap-annot \
    --print-coefficients \
    --print-delete-vals \
    --out ${OUT_DIR}/${TRAIT}_partitioned

# Output ${TRAIT}_partitioned.results columns:
#   Category | Prop._SNPs | Prop._h2 | Prop._h2_std_error | Enrichment | Enrichment_std_error | Enrichment_p
# Apply Bonferroni: 0.05 / nrow(partitioned.results) (~ 0.05/97 = 5.2e-4 for baseline-LD v2.2)

# ============================================================
# Step 4: Cell-type / tissue prioritization (Finucane 2018)
# ============================================================
# Identifies trait-relevant tissues by testing per-cell-type chromatin annotations
# marginal to the baseline model.

ldsc.py \
    --h2-cts ${MUNGED} \
    --ref-ld-chr ${REF_DIR}/1000G_Phase3_baseline/baseline. \
    --ref-ld-chr-cts ${CTS_MANIFEST} \
    --w-ld-chr ${WEIGHTS} \
    --out ${OUT_DIR}/${TRAIT}_cts

# Output ${TRAIT}_cts.cell_type_results.txt:
#   Name | Coefficient | Coefficient_std_error | Coefficient_P_value
# Apply Bonferroni: 0.05 / N_tissues (~ 2.5e-4 for the ~200-tissue manifest)
# Top trait-relevant tissues: positive Coefficient AND p < threshold

# ============================================================
# Step 5: Cross-trait genetic correlation (rg) -- example with a paired sumstats file
# ============================================================
# Uncomment and supply a second munged sumstats file to run rg
#
# SECOND_MUNGED=${OUT_DIR}/${TRAIT2}.sumstats.gz
# ldsc.py \
#     --rg ${MUNGED},${SECOND_MUNGED} \
#     --ref-ld-chr ${EUR_LD}/ \
#     --w-ld-chr ${EUR_LD}/ \
#     --out ${OUT_DIR}/${TRAIT}_vs_${TRAIT2}_rg
#
# Output reports: rg, se, p, h2_obs, h2_int, gcov_int
# gcov_int absorbs sample overlap; rg is unbiased even with overlap

echo "LDSC pipeline complete. Results in ${OUT_DIR}/"
echo "  ${TRAIT}_h2.log               -> total h2 + intercept diagnostics"
echo "  ${TRAIT}_partitioned.results  -> baseline-LD v2.2 enrichment"
echo "  ${TRAIT}_cts.cell_type_results.txt -> tissue prioritization"

#!/bin/bash
# Reference: LDAK 6.0+, LDAK-Thin + BaselineLD tagging files | Verify API if version differs
#
# LDAK SumHer alternative h2 / enrichment estimate for reconciliation with LDSC.
# Speed 2019 Nat Genet 51:277 introduces the LDAK-Thin model (MAF + LD reweighting).
# Gazal 2019 Nat Genet 51:1202 reconciles LDSC vs LDAK enrichment differences.
#
# Operational rule: when functional enrichment is the primary claim, report BOTH
# LDSC baseline-LD and LDAK SumHer; treat > 2x discordance as model-dependent.
#
# Usage:
#   bash ldak_sumher.sh <gwas_ldak.txt> <trait_prefix>
#
# Input format (LDAK requires header):
#   Predictor A1 A2 n Z
#   rs123 A G 100000 -2.34
#   ...
# Where Predictor = SNP ID, A1 = effect allele, A2 = other allele,
#       n = per-SNP sample size, Z = signed Z-score.
# Convert from BETA/SE: Z = BETA / SE.

set -euo pipefail

GWAS_LDAK=${1:?Usage: ldak_sumher.sh <gwas_ldak.txt> <trait_prefix>}
TRAIT=${2:?provide trait_prefix}

LDAK=${LDAK_BIN:-./ldak6.3.linux}
TAGFILE=${LDAK_TAGFILE:-./ldak.thin.hapmap.gbr.tagging}
# Pre-computed LDAK-Thin tagging files at https://dougspeed.com/pre-computed-tagging-files/
# ldak.thin.hapmap.gbr.tagging covers HapMap3 + GBR ancestry; matches LDSC's HapMap3 default

OUT_DIR=${TRAIT}_ldak
mkdir -p ${OUT_DIR}

# ============================================================
# Step 1: Total h2 estimate under LDAK-Thin model
# ============================================================
# --check-sums NO skips strict sumstat consistency checks (use only if GWAS pre-QC'd)

${LDAK} --sum-hers ${OUT_DIR}/${TRAIT}_sumher \
    --summary ${GWAS_LDAK} \
    --tagfile ${TAGFILE} \
    --check-sums NO

# Output ${TRAIT}_sumher.hers:
#   Component | Heritability | Std_Error | Influence | Z-Score | P-Value
# The single-component output is total h2 under LDAK-Thin

# ============================================================
# Step 2: Partitioned h2 with BaselineLD annotations
# ============================================================
# The BaselineLD model used here supplies 86 LD/MAF/functional annotations (binary + continuous) via --annotation-number 86.
# BLD-LDAK is a distinct standalone 65-annotation model (LD/MAF bins), NOT additive on top of BaselineLD.
# Pre-computed BaselineLD annotations at dougspeed.com/resources (download BaselineLD.zip =
# 96 annotations; extract to ./BaselineLD/BaselineLD1 ... BaselineLD96, the run uses the first 86).
# The annotation flags belong to --calc-tagging (which builds the tagging file from a
# genotype reference), NOT to --sum-hers. LDAK uses --annotation-number/--annotation-prefix
# (continuous annotations) or --partition-number/--partition-prefix (binary partitions) --
# there is no --catfile flag. Partitioned h2 is therefore a two-step workflow: build an
# annotated tagging file, then estimate h2 against it with no annotation flags on --sum-hers.

BLD_ANNOT_PREFIX=${LDAK_ANNOT_PREFIX:-./BaselineLD/BaselineLD}
N_BLD_ANNOTS=${LDAK_N_ANNOTS:-86}
REF_BFILE=${LDAK_REF_BFILE:-./ref_panel}

# 2a. Build the annotated tagging file (--power -.25 sets the assumed heritability model)
${LDAK} --calc-tagging ${OUT_DIR}/${TRAIT}_bld \
    --bfile ${REF_BFILE} \
    --power -.25 \
    --annotation-number ${N_BLD_ANNOTS} \
    --annotation-prefix ${BLD_ANNOT_PREFIX}

# 2b. Estimate partitioned h2 against that tagging file (no annotation flags here)
${LDAK} --sum-hers ${OUT_DIR}/${TRAIT}_bld \
    --summary ${GWAS_LDAK} \
    --tagfile ${OUT_DIR}/${TRAIT}_bld.tagging \
    --check-sums NO

# Output files:
#   ${TRAIT}_bld.hers    -> per-component heritability + SE
#   ${TRAIT}_bld.enrich  -> per-category enrichment + Z-Score + P-Value
# Apply Bonferroni at 0.05 / N_categories for per-category claims

# ============================================================
# Step 3: LDSC vs LDAK reconciliation report
# ============================================================
# Cross-reference with the LDSC partitioned output (see ldsc_partitioned_h2.sh).
# Per Gazal 2019 Nat Genet 51:1202:
#   - For functional enrichment claims, report BOTH model estimates
#   - If LDSC and LDAK disagree by > 2x, flag as model-dependent
#   - Prefer LDAK SumHer for conserved-region enrichment per Speed 2019

echo "LDAK SumHer pipeline complete. Results in ${OUT_DIR}/"
echo "  ${TRAIT}_sumher.hers   -> total h2 (LDAK-Thin)"
echo "  ${TRAIT}_bld.hers      -> per-component h2 contributions"
echo "  ${TRAIT}_bld.enrich    -> BaselineLD per-category enrichment + Z + P"
echo ""
echo "Compare against LDSC partitioned output:"
echo "  - Total h2: LDSC GCTA model vs LDAK-Thin (typically within 20%)"
echo "  - Per-category enrichment: discordance > 2x indicates model-dependence"
echo "  - Cite Gazal 2019 Nat Genet 51:1202 in methods when reporting both"

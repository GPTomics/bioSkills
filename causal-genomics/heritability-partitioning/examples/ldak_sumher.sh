#!/bin/bash
# Reference: LDAK 6.0+, BLD-LDAK tagging files | Verify API if version differs
#
# LDAK SumHer alternative h2 / enrichment estimate for reconciliation with LDSC.
# Speed 2019 Nat Genet 51:277 introduces the LDAK-Thin model (MAF + LD reweighting).
# Hou 2019 Nat Genet 51:1244 reconciles LDSC vs LDAK enrichment differences.
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

LDAK=${LDAK_BIN:-./ldak6.linux}
TAGFILE=${LDAK_TAGFILE:-./bld.ldak.thin.hapmap.gbr.tagging}
# Pre-computed BLD-LDAK tagging files at https://dougspeed.com/pre-computed-tagging-files/
# bld.ldak.thin.hapmap.gbr.tagging covers HapMap3 + GBR ancestry; matches LDSC's HapMap3 default

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
# Step 2: Partitioned h2 with BLD-LDAK annotations
# ============================================================
# BLD-LDAK adds 65 LD- and MAF-stratified bins on top of the BaselineLD annotations.
# Pre-computed BaselineLD annotations available at dougspeed.com/bldldak (download
# BaselineLD.zip and extract to ./BaselineLD/BaselineLD1 ... BaselineLD86).
# LDAK uses --annotation-number/--annotation-prefix (continuous annotations) or
# --partition-number/--partition-prefix (binary partitions) -- there is no --catfile flag.
# When the tagging file already encodes the BLD-LDAK categories (e.g.
# bld.ldak.thin.hapmap.gbr.tagging), --sum-hers automatically reports per-category
# enrichment in the output .enrich file without extra annotation flags.

BLD_ANNOT_PREFIX=${LDAK_ANNOT_PREFIX:-./BaselineLD/BaselineLD}
N_BLD_ANNOTS=${LDAK_N_ANNOTS:-86}

${LDAK} --sum-hers ${OUT_DIR}/${TRAIT}_bld \
    --summary ${GWAS_LDAK} \
    --tagfile ${TAGFILE} \
    --annotation-number ${N_BLD_ANNOTS} \
    --annotation-prefix ${BLD_ANNOT_PREFIX} \
    --check-sums NO

# Output files:
#   ${TRAIT}_bld.hers    -> per-component heritability + SE
#   ${TRAIT}_bld.enrich  -> per-category enrichment + Z-Score + P-Value
# Apply Bonferroni at 0.05 / N_categories for per-category claims

# ============================================================
# Step 3: LDSC vs LDAK reconciliation report
# ============================================================
# Cross-reference with the LDSC partitioned output (see ldsc_partitioned_h2.sh).
# Per Hou 2019 Nat Genet 51:1244:
#   - For functional enrichment claims, report BOTH model estimates
#   - If LDSC and LDAK disagree by > 2x, flag as model-dependent
#   - Prefer LDAK SumHer for conserved-region enrichment per Speed 2019

echo "LDAK SumHer pipeline complete. Results in ${OUT_DIR}/"
echo "  ${TRAIT}_sumher.hers   -> total h2 (LDAK-Thin)"
echo "  ${TRAIT}_bld.hers      -> per-component h2 contributions"
echo "  ${TRAIT}_bld.enrich    -> BLD-LDAK per-category enrichment + Z + P"
echo ""
echo "Compare against LDSC partitioned output:"
echo "  - Total h2: LDSC GCTA model vs LDAK-Thin (typically within 20%)"
echo "  - Per-category enrichment: discordance > 2x indicates model-dependence"
echo "  - Cite Hou 2019 Nat Genet 51:1244 in methods when reporting both"

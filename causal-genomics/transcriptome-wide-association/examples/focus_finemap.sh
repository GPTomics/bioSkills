#!/usr/bin/env bash
# Reference: MetaXcan 0.7+, FUSION 2.0+, FOCUS 0.8+, plink2 | Verify API if version differs
# FOCUS probabilistic gene-level fine-mapping of TWAS hits.
# Resolves co-significant gene clusters at gene-dense loci into a credible causal-gene set
# with per-gene posterior inclusion probabilities (PIPs).

set -euo pipefail

# ---- Inputs ----
GWAS_FILE='gwas.sumstats'                                # GWAS sumstats (CHR SNP BP A1 A2 Z P columns at minimum)
LD_REF_PREFIX='1000G_EUR/chr'                            # PLINK bfile per chromosome (chr1.bim/bed/fam, ...)
FOCUS_DB='focus_gtex_v8_whole_blood.db'                  # FOCUS DB matched to the TWAS weight panel
TISSUE='Whole_Blood'                                     # Tissue label inside the FOCUS DB

# Genome-wide significance threshold for SNPs that flag a locus for fine-mapping.
# 5e-8 is standard GWAS genome-wide significance; loci below this threshold trigger FOCUS.
P_THRESHOLD='5e-8'

# Locus window: FOCUS extracts SNPs within +/- this distance of each genome-wide
# significant lead. 1 Mb captures essentially all cis-eQTL signal and downstream LD.
LOCUS_WINDOW='1000000'

OUT_PREFIX='gwas_focus_whole_blood'

# ---- Step 1: build FOCUS database from FUSION weights if not pre-built ----
# Skip this if using a pre-built FOCUS DB. Custom panels need this step.
# focus import-fusion gtex_whole_blood.pos --output focus_gtex_v8_whole_blood.db

# ---- Step 2: run FOCUS fine-mapping ----
# FOCUS takes the underlying GWAS sumstats (not the TWAS Z output) and internally:
#   1. computes per-gene TWAS Z using the panel weights
#   2. estimates the gene-by-gene predicted-expression correlation matrix as the LD analog
#   3. runs a Bayesian variable-selection model over genes (analog of variant fine-mapping)
focus finemap \
    "${GWAS_FILE}" \
    "${LD_REF_PREFIX}" \
    "${FOCUS_DB}" \
    --p-threshold "${P_THRESHOLD}" \
    --tissue "${TISSUE}" \
    --locus-window "${LOCUS_WINDOW}" \
    --out "${OUT_PREFIX}"
# Single-ancestry mode does not take --locations; that flag is MA-FOCUS-only
# (signaled by colon-separated per-ancestry inputs + --locations 38:EUR-EAS-AFR).

# ---- Step 3: filter credible-set genes ----
# FOCUS reports per-gene PIP. PIP >= 0.8 = causal candidate; 0.5 <= PIP < 0.8 = suggestive;
# PIP < 0.5 at a co-significant locus = LD-tagged co-regulated gene (NOT causal).
# Mancuso 2019 Nat Genet 51:675 convention.
PIP_THRESH='0.8'

awk -v t="${PIP_THRESH}" -F'\t' '
    NR==1 {print; for (i=1; i<=NF; i++) col[$i]=i; next}
    $col["pip"] >= t {print}
' "${OUT_PREFIX}.focus.tsv" > "${OUT_PREFIX}.credible.tsv"

echo 'FOCUS fine-mapping complete.'
echo "Full output: ${OUT_PREFIX}.focus.tsv"
echo "Causal-gene candidates (PIP >= ${PIP_THRESH}): ${OUT_PREFIX}.credible.tsv"

# ---- Step 4 (optional): MA-FOCUS for multi-ancestry ----
# Requires per-ancestry sumstats, per-ancestry LD references, and per-ancestry FOCUS DBs.
# Joint inference assumes a shared causal gene across ancestries; trans-ethnic
# gene-effect heterogeneity violates this and produces inflated heterogeneous-group probability.
#
# MA-FOCUS reuses the single-ancestry `focus finemap` CLI; multi-ancestry mode
# is signaled by colon-separated per-ancestry sumstats/LD/weight DBs plus
# paired ancestry codes in --locations.
# focus finemap \
#     eur.sumstats.tsv.gz:eas.sumstats.tsv.gz:afr.sumstats.tsv.gz \
#     1000G_EUR/chr:1000G_EAS/chr:1000G_AFR/chr \
#     focus_eur.db:focus_eas.db:focus_afr.db \
#     --p-threshold "${P_THRESHOLD}" \
#     --tissue "${TISSUE}" \
#     --locations 38:EUR-EAS-AFR \
#     --out gwas_ma_focus

#!/usr/bin/env bash
# Reference: SuSiEx head (getian107/SuSiEx 2024), PLINK 1.9+ | Verify CLI flags with `SuSiEx --help` if version differs
##
## Cross-ancestry joint fine-mapping with SuSiEx.
## Assumes per-population GWAS summary statistics and reference panels are pre-built.
## Output: shared credible sets that exploit population-specific LD to shrink set size.

set -euo pipefail

LOCUS_NAME=locus1
LOCUS_CHR=6
LOCUS_BP_START=30000000
LOCUS_BP_END=31000000

EUR_SST=eur_sumstats.txt
EAS_SST=eas_sumstats.txt
AFR_SST=afr_sumstats.txt

EUR_REF=1000G_EUR
EAS_REF=1000G_EAS
AFR_REF=1000G_AFR

N_EUR=500000
N_EAS=200000
N_AFR=80000

EUR_LD=eur_ld
EAS_LD=eas_ld
AFR_LD=afr_ld

OUT_DIR=susiex_out
mkdir -p "${OUT_DIR}"

# Each summary-stat file requires columns CHR/SNP/BP/A1/A2/BETA/SE/P (column numbers
# supplied via --chr_col, --snp_col, --bp_col, --a1_col, --a2_col, --eff_col, --se_col,
# --pval_col). Each reference panel: PLINK bim/bed/fam triplet covering the locus.
# --ld_file is required even when --ref_file is provided; LD matrices can be pre-built
# with `plink --r square` per population, or use SuSiEx's bundled `getLD` helper.

SuSiEx \
    --sst_file="${EUR_SST},${EAS_SST},${AFR_SST}" \
    --n_gwas="${N_EUR},${N_EAS},${N_AFR}" \
    --ref_file="${EUR_REF},${EAS_REF},${AFR_REF}" \
    --ld_file="${EUR_LD},${EAS_LD},${AFR_LD}" \
    --chr="${LOCUS_CHR}" \
    --bp="${LOCUS_BP_START},${LOCUS_BP_END}" \
    --chr_col=1,1,1 \
    --snp_col=2,2,2 \
    --bp_col=3,3,3 \
    --a1_col=4,4,4 \
    --a2_col=5,5,5 \
    --eff_col=6,6,6 \
    --se_col=7,7,7 \
    --pval_col=8,8,8 \
    --pop=EUR,EAS,AFR \
    --out_dir="${OUT_DIR}" \
    --out_name="${LOCUS_NAME}" \
    --level=0.95 \
    --threads=4

# Outputs:
#   ${OUT_DIR}/${LOCUS_NAME}.snp   per-SNP PIP and credible-set membership (joint)
#   ${OUT_DIR}/${LOCUS_NAME}.cs    credible set summary (size, purity, top SNP)
#   ${OUT_DIR}/${LOCUS_NAME}.log   convergence diagnostics

# Inspect the credible set summary:
cat "${OUT_DIR}/${LOCUS_NAME}.cs"

# Compare to single-ancestry susie_rss (EUR-only): typical observation is that
# credible-set size shrinks 2-5x when AFR (shorter LD blocks) is included.
# Document the per-population sample sizes; SuSiEx assumes shared causal variants.

echo
echo "Reference: Yuan K et al 2024 Nat Genet 56:1841 (SuSiEx)"

#!/usr/bin/env bash
# Reference: MetaXcan 0.7+, FUSION 2.0+, FOCUS 0.8+, plink2 | Verify API if version differs
# S-PrediXcan + S-MultiXcan pipeline: run TWAS across all GTEx v8 MASHR-EUR tissues
# and combine via the joint multi-tissue test.

set -euo pipefail

# ---- Inputs ----
GWAS_FILE='gwas.txt'                                              # harmonised GWAS sumstats
MODELS_DIR='mashr_models'                                         # PredictDB GTEx v8 MASHR-EUR (.db + .txt.gz per tissue)
SNP_COV='gtex_v8_expression_mashr_snp_covariance.txt.gz'          # MASHR cross-tissue SNP covariance
OUT_DIR='spredixcan_out'

mkdir -p "${OUT_DIR}"

# ---- Tissues to run ----
# GTEx v8 has 49 tissues; this list is representative. Edit as needed.
# Tissues with N < 100 donors are flagged as unstable (e.g. some brain sub-regions); skip
# them for biology-agnostic screens.
TISSUES=(
    Whole_Blood
    Liver
    Adipose_Subcutaneous
    Adipose_Visceral_Omentum
    Brain_Frontal_Cortex_BA9
    Brain_Cortex
    Brain_Hippocampus
    Brain_Cerebellum
    Heart_Left_Ventricle
    Heart_Atrial_Appendage
    Muscle_Skeletal
    Pancreas
    Lung
    Thyroid
)

# ---- Step 1: per-tissue S-PrediXcan ----
for tissue in "${TISSUES[@]}"; do
    echo "Running S-PrediXcan: ${tissue}"
    python SPrediXcan.py \
        --model_db_path "${MODELS_DIR}/mashr_${tissue}.db" \
        --covariance "${MODELS_DIR}/mashr_${tissue}.txt.gz" \
        --gwas_file "${GWAS_FILE}" \
        --snp_column SNP \
        --effect_allele_column A1 \
        --non_effect_allele_column A2 \
        --beta_column BETA \
        --pvalue_column P \
        --output_file "${OUT_DIR}/${tissue}.csv"
done

# ---- Step 2: S-MultiXcan joint multi-tissue ----
# Joint test combines per-tissue Z-scores via PCA-regularised regression on the inter-tissue
# correlation matrix. Regularisation 0.1 (default) is conservative; raise to 0.5 if tissues
# are near-collinear (e.g. multiple brain sub-regions).
python SMulTiXcan.py \
    --models_folder "${MODELS_DIR}" \
    --models_name_pattern 'mashr_(.*)\.db' \
    --snp_covariance "${SNP_COV}" \
    --metaxcan_folder "${OUT_DIR}" \
    --metaxcan_filter '(.*)\.csv' \
    --metaxcan_file_name_parse_pattern '(.*)\.csv' \
    --gwas_file "${GWAS_FILE}" \
    --snp_column SNP \
    --effect_allele_column A1 \
    --non_effect_allele_column A2 \
    --beta_column BETA \
    --pvalue_column P \
    --regularization 0.1 \
    --cutoff_condition_number 30 \
    --verbosity 9 \
    --throw \
    --output 'joint_multitissue.csv'

# ---- Step 3: filter significant genes ----
# Genome-wide TWAS Bonferroni: 0.05 / 22000 ~ 2.3e-6
# S-MultiXcan returns one p per gene (across all tissues jointly); pvalue is column 'pvalue'.
# Resolve its column index from the header instead of assuming $NF, since SMulTiXcan also
# emits z_min/z_max/eigen_*/tmi/status columns after pvalue.
awk -F',' 'NR==1 { for (i=1; i<=NF; i++) if ($i == "pvalue") pcol = i; print; next }
           pcol && $pcol != "" && $pcol+0 < 2.3e-6' joint_multitissue.csv > joint_multitissue_sig.csv

echo 'S-PrediXcan + S-MultiXcan pipeline complete.'
echo "Per-tissue outputs: ${OUT_DIR}/<tissue>.csv"
echo 'Joint multi-tissue: joint_multitissue.csv'
echo 'Significant genes (p < 2.3e-6): joint_multitissue_sig.csv'

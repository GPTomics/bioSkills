#!/usr/bin/env bash
# Reference: chemprop 2.0+, RDKit 2024.09+ | Verify API if version differs
# QSAR pipeline: scaffold-split, 5-fold x 5 ensemble chemprop 2.x D-MPNN + RDKit 2D features
# chemprop 2.x: 'chemprop train' (space; dashes not underscores). 1.x legacy: 'chemprop_train ...'

set -euo pipefail

DATA_PATH="${1:-data.csv}"
SAVE_DIR="${2:-chemprop_model}"
TASK_TYPE="${3:-classification}"

# scaffold_balanced split prevents chemotype leakage; mandatory for production QSAR
# Metric 'roc' for binary classification; 'mae' / 'rmse' for regression
chemprop train \
    --data-path "${DATA_PATH}" \
    --task-type "${TASK_TYPE}" \
    --save-dir "${SAVE_DIR}" \
    --molecule-featurizers rdkit_2d_normalized \
    --num-folds 5 \
    --ensemble-size 5 \
    --epochs 50 \
    --batch-size 128 \
    --split scaffold_balanced \
    --split-sizes 0.8 0.1 0.1 \
    --metric roc \
    --data-seed 42

# Predict on new SMILES with ensemble uncertainty
chemprop predict \
    --test-path new_compounds.csv \
    --model-path "${SAVE_DIR}/best.pt" \
    --preds-path predictions.csv \
    --molecule-featurizers rdkit_2d_normalized

# Output: predictions.csv with columns <target>_pred + <target>_unc (ensemble variance)
# High variance compounds are out-of-distribution; predictions are less reliable

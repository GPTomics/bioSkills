#!/usr/bin/env bash
# Reference: chemprop 2.2.x, RDKit 2024.09+ | Verify API if version differs
# QSAR pipeline: five scaffold-split replicates, each with a five-model chemprop 2.x ensemble
# chemprop 2.x: 'chemprop train' (space; dashes not underscores). 1.x legacy: 'chemprop_train ...'

set -euo pipefail

DATA_PATH="${1:-data.csv}"
SAVE_DIR="${2:-chemprop_model}"
TASK_TYPE="${3:-classification}"

case "${TASK_TYPE}" in
    classification)
        METRIC="roc"
        ;;
    regression)
        METRIC="rmse"
        ;;
    *)
        echo "Unsupported task type: ${TASK_TYPE}" >&2
        exit 2
        ;;
esac

# scaffold_balanced keeps scaffold groups together; choose it when that matches deployment
# Metric 'roc' for binary classification; 'mae' / 'rmse' for regression
chemprop train \
    --data-path "${DATA_PATH}" \
    --task-type "${TASK_TYPE}" \
    --save-dir "${SAVE_DIR}" \
    --molecule-featurizers rdkit_2d \
    --num-replicates 5 \
    --ensemble-size 5 \
    --epochs 50 \
    --batch-size 128 \
    --split scaffold_balanced \
    --split-sizes 0.8 0.1 0.1 \
    --metric "${METRIC}" \
    --data-seed 42

# Prediction requires the actual checkpoint path(s) produced by the installed
# chemprop version. Inspect `chemprop predict --help` and the saved run directory;
# do not assume a `${SAVE_DIR}/best.pt` layout or automatic `_unc` columns.
# To request ensemble uncertainty, supply the emitted checkpoint paths and add
# `--uncertainty-method ensemble` to `chemprop predict`.

# Classification Models for Omics Data Usage Guide

## Overview

Build diagnostic and prognostic classifiers on omics feature matrices with regularized logistic regression, random forest, and gradient-boosted trees. The discipline this skill enforces is that in the p>>n regime a regularized linear model is often the strongest baseline, the probability (not the label) is usually the product so calibration matters more than accuracy, and a suspiciously perfect AUC is most often a batch artifact rather than biology.

## Prerequisites

```bash
pip install scikit-learn xgboost imbalanced-learn pandas
```

Conceptual prerequisites: know whether the deliverable is a probability or a hard label, whether batch is confounded with the outcome, and whether feature selection and scaling are inside the cross-validation loop (machine-learning/model-validation).

## Quick Start

Tell your AI agent what you want to do:
- "Build an elastic-net logistic classifier as a baseline on my expression data"
- "Train XGBoost with early stopping and check whether it beats the linear model"
- "Is my high AUC real, or is it learning the batch?"
- "My classes are imbalanced -- should I use SMOTE?"

## Example Prompts

### Algorithm choice

> "Build a classifier from my RNA-seq matrix. Start with regularized logistic regression and only escalate to trees if there is interaction signal it misses."

> "Train an XGBoost classifier with a low learning rate, shallow depth, and early stopping; the dataset is small."

### Batch shortcuts

> "My cross-validated AUC is 0.97. Check whether the model is learning batch instead of biology by predicting the batch from the features and using a batch-aware split."

### Imbalance and calibration

> "My positive class is 8%. Handle the imbalance without destroying calibration, and tell me whether SMOTE is appropriate."

> "Check whether my random forest probabilities are calibrated and recalibrate them if not."

### Mixed features and missing data

> "My features mix continuous expression and categorical genotype with missing values. Set up the preprocessing and a model that handles NaN natively."

## What the Agent Will Do

1. Start with a regularized linear baseline (often the ceiling in p>>n)
2. Escalate to RF/XGBoost only when nonlinear or interaction signal is present, with early stopping
3. Run the batch-shortcut checks (predict the batch; leave-one-batch-out) before trusting any AUC
4. Handle imbalance by class weighting or threshold tuning, not resampling, for risk models
5. Check and recalibrate probabilities; report AUPRC and MCC under imbalance
6. Keep all preprocessing inside the CV fold and defer unbiased evaluation to model-validation

## Tips

- "Random forest is the obvious choice for expression data" is a myth; regularized logistic/linear SVM frequently win in p>>n (Statnikov 2008)
- GBDTs beat deep nets on tabular/omics matrices (Grinsztajn 2022); reserve deep nets for very large n or raw/multimodal inputs
- RF compresses probabilities toward 0.5, boosting pushes them to the extremes -- both need recalibration if the probability is used
- SMOTE/oversampling destroys calibration for no AUC gain on risk models; class-weight or tune the threshold instead, and if you must resample use an `imblearn` Pipeline
- XGBoost 2.x: set `early_stopping_rounds` and `eval_metric` in the constructor, not `fit()`
- XGBoost/LightGBM handle missing values natively; in omics, missingness is often informative (below detection limit)
- Trees are scale-invariant; standardize for logistic/SVM/kNN/PCA

## Related Skills

- machine-learning/model-validation - Nested CV, calibration, and net benefit for the trained classifier
- machine-learning/biomarker-discovery - Select features before modeling (inside the CV fold)
- machine-learning/prediction-explanation - Interpret the classifier and detect shortcuts with SHAP
- machine-learning/survival-analysis - Time-to-event outcomes that classifiers cannot handle
- differential-expression/batch-correction - Batch correction done design-aware, not across the split
- expression-matrix/normalization - Per-sample normalization that is safe outside the CV fold

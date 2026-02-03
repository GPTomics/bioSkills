# Classification Models Usage Guide

## Overview

Build classification models for biomarker discovery and diagnostics using RandomForest, XGBoost, and logistic regression with sklearn-compatible APIs.

## Prerequisites

```bash
pip install scikit-learn xgboost pandas matplotlib
```

## Quick Start

Tell your AI agent what you want to do:
- "Build a random forest classifier on my expression data"
- "Train an XGBoost model to predict disease status from my omics data"
- "Use logistic regression with L1 penalty to find sparse biomarkers"
- "Evaluate my classifier with ROC curves and AUC"

## Example Prompts

### Basic Classification

> "Train a random forest classifier on my gene expression matrix (expression.csv) with labels in metadata.csv. Show me the ROC curve and top 20 important features."

> "Build an XGBoost classifier to predict tumor vs normal from my RNA-seq counts."

### Biomarker Selection

> "Use L1-regularized logistic regression to find a sparse set of genes that classify my samples. Which genes have non-zero coefficients?"

### Multi-class

> "I have expression data with 4 cancer subtypes. Train a classifier and show the confusion matrix."

### Model Comparison

> "Compare RandomForest, XGBoost, and logistic regression on my dataset. Which has the best AUC?"

## What the Agent Will Do

1. Load expression matrix and sample labels
2. Split data into train/test sets (stratified)
3. Build preprocessing pipeline with scaling
4. Train classifier and generate predictions
5. Report metrics (accuracy, AUC, classification report)
6. Visualize ROC curve and feature importances

## Tips

- Always use stratified splits for imbalanced classes
- StandardScaler is essential for logistic regression; optional for tree models
- Use class_weight='balanced' when classes are imbalanced
- XGBoost parameters: use random_state not seed, n_jobs not nthread
- L1 penalty gives sparse solutions; L2 handles correlated features better
- For small datasets, consider nested cross-validation instead of train/test split

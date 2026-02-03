# Biomarker Pipeline Usage Guide

## Overview

End-to-end workflow for biomarker discovery combining feature selection, model training with nested cross-validation, interpretation, and validation. Produces a validated biomarker panel with an accompanying classifier.

## Prerequisites

```bash
pip install scikit-learn boruta shap xgboost pandas numpy matplotlib joblib
```

**Input data:**
- Expression matrix (genes x samples) as CSV
- Metadata file with sample IDs and condition/label column

## Quick Start

Tell your AI agent what you want to do:
- "Build a biomarker classifier from my expression data"
- "Run the full biomarker discovery pipeline with nested CV"
- "Select features and train a validated classifier"
- "Create a biomarker panel for disease vs control classification"

## Example Prompts

### Basic Biomarker Discovery

> "I have expression.csv and metadata.csv. Build a biomarker classifier for my disease vs control samples."

> "Run biomarker discovery with LASSO stability selection and nested cross-validation."

### With Specific Methods

> "Use Boruta for feature selection and train a Random Forest classifier with SHAP interpretation."

> "Build a minimal biomarker signature using LASSO with strict stability threshold (0.8)."

### Validation Focus

> "Create a validated biomarker panel with bootstrap confidence intervals for AUC."

> "Train a classifier with nested CV and export the model for external validation."

## What the Agent Will Do

1. Load and prepare data with stratified train/test split
2. Scale features (fit on training only to prevent leakage)
3. Select features using Boruta or LASSO with stability selection
4. Train classifier with nested CV for unbiased performance estimation
5. Generate SHAP plots for model interpretation
6. Validate on held-out test set with bootstrap confidence intervals
7. Export biomarker panel and trained model

## Tips

- Start with at least 20 samples per class for reasonable statistical power
- Use Boruta for comprehensive biomarker panels (finds all relevant features)
- Use LASSO for minimal signatures (finds sparse feature sets)
- Always use nested CV to avoid overfitting bias in performance estimates
- Check that SHAP top features align with selected features (sanity check)
- Pre-filter with differential expression if starting with >10k features
- Consider biological validation with independent dataset or orthogonal assay
- Class imbalance >3:1 may require stratification or SMOTE

## Related Skills

- machine-learning/biomarker-discovery - Detailed feature selection methods
- machine-learning/model-validation - Nested CV implementation details
- machine-learning/omics-classifiers - Classifier options and tuning
- machine-learning/prediction-explanation - SHAP and LIME interpretation
- differential-expression/de-results - Pre-filter with DE genes
- pathway-analysis/go-enrichment - Functional enrichment of biomarkers

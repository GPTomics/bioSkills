# Biomarker Pipeline Usage Guide

## Overview

End-to-end workflow for biomarker discovery combining feature selection, leakage-safe cross-validation, calibration, interpretation, and validation. Produces an honestly-validated biomarker panel with an accompanying classifier.

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
- "Run the full biomarker discovery pipeline with leakage-safe cross-validation"
- "Select features and train a validated classifier"
- "Create a biomarker panel for disease vs control classification"

## Example Prompts

### Basic Biomarker Discovery

> "I have expression.csv and metadata.csv. Build a biomarker classifier for my disease vs control samples."

> "Run biomarker discovery with LASSO stability selection and leakage-safe cross-validation."

### With Specific Methods

> "Use Boruta for feature selection and train a Random Forest classifier with SHAP interpretation."

> "Build a minimal biomarker signature using LASSO with strict stability threshold (0.8)."

### Validation Focus

> "Create a validated biomarker panel with bootstrap confidence intervals for AUC."

> "Train a classifier with leakage-safe cross-validation and export the model for external validation."

### Pipeline-level decisions

> "I have multiple biopsies/timepoints per patient - how do I split so I don't leak subjects?"

> "My AUC is great in CV but collapses on an external cohort - what leaked?"

> "The panel will output risk scores - do I need calibration on top of AUC?"

## What the Agent Will Do

1. Load and prepare data with a group- and class-aware split by subject (no subject in both train and test)
2. Scale features (fit on training only to prevent leakage)
3. Select features using Boruta or LASSO with stability selection
4. Estimate performance with leakage-safe cross-validation (selection inside each fold)
5. Audit the model with interventional SHAP (shortcut/batch detection), aggregated over modules
6. Validate on held-out test set with bootstrap confidence intervals and calibration
7. Export biomarker panel and trained model

## Tips

- Start with at least 20 samples per class for reasonable statistical power
- Use Boruta for comprehensive biomarker panels (finds all relevant features)
- Use LASSO for minimal signatures (finds sparse feature sets)
- Split by the SUBJECT (patient/donor/site), not the sample; repeated biopsies/timepoints from one subject in both train and test is group leakage a held-out set cannot detect
- Keep feature selection inside the CV fold; selection-before-CV inflates AUC toward 1.0 even on noise
- SHAP is an audit for batch/housekeeping shortcuts, not a check that it matches selection -- divergence under correlated features is expected
- Pre-filter with differential expression if starting with >10k features
- Consider biological validation with independent dataset or orthogonal assay
- For class imbalance prefer `class_weight='balanced'` and threshold-moving; do not SMOTE risk models (it inflates minority probabilities and breaks calibration)

## Related Skills

- database-access/geo-data - Public expression cohorts for independent validation
- database-access/sra-data - Pull raw FASTQ to build re-quantified validation cohorts
- database-access/uniprot-access - Protein-level features (sequence, GO, PTMs) for protein biomarkers
- machine-learning/biomarker-discovery - Detailed feature selection methods
- machine-learning/model-validation - Nested CV implementation details
- machine-learning/omics-classifiers - Classifier options and tuning
- machine-learning/prediction-explanation - SHAP and LIME interpretation
- differential-expression/de-results - Pre-filter with DE genes
- pathway-analysis/go-enrichment - Functional enrichment of biomarkers

# Cross-Validation Usage Guide

## Overview

Implement proper cross-validation strategies to get unbiased performance estimates on biomedical datasets and avoid overfitting during biomarker discovery.

## Prerequisites

```bash
pip install scikit-learn pandas numpy
```

## Quick Start

Tell your AI agent what you want to do:
- "Run nested cross-validation on my classifier"
- "Evaluate my model with stratified 5-fold CV"
- "Use leave-one-out validation for my small dataset"
- "Make sure patients from the same individual are kept together in CV"

## Example Prompts

### Nested CV

> "Run nested cross-validation on my expression classifier. Use 5 outer folds for evaluation and 3 inner folds for hyperparameter tuning with grid search."

> "I want an unbiased estimate of my classifier's performance. Run nested CV and report the mean and standard deviation of AUC."

### Stratified CV

> "Evaluate my random forest with stratified 5-fold cross-validation. My classes are imbalanced."

### Small Datasets

> "I only have 25 samples. Run leave-one-out cross-validation to maximize training data."

### Grouped Data

> "My dataset has multiple samples per patient. Run cross-validation keeping all samples from each patient in the same fold."

## What the Agent Will Do

1. Identify appropriate CV strategy for dataset size
2. Set up outer CV for evaluation (and inner CV if nested)
3. Ensure feature selection/preprocessing is inside the CV loop
4. Calculate performance metrics per fold
5. Report mean +/- std of performance

## Tips

- Always use stratified splits when classes are imbalanced
- Put ALL preprocessing (scaling, feature selection) inside the CV loop to avoid data leakage
- Nested CV is essential when tuning hyperparameters on small datasets
- Use GroupKFold when samples aren't independent (repeated measures, same patient)
- Leave-one-out gives maximum training data but high variance estimates
- Report standard deviation of CV scores, not just the mean

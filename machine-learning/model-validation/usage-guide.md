# Model Validation Usage Guide

## Overview

Estimate the performance of predictive models on omics and biomedical data honestly: keep every data-dependent step inside the resampling loop, match the cross-validation scheme to the real generalization question, and report calibration and net benefit rather than discrimination alone. The two failure roots this skill guards against are leakage (test information entering model construction) and grading a model on the same data used to choose it.

## Prerequisites

```bash
pip install scikit-learn pandas numpy
```

Conceptual prerequisites: know the unit of independence (patient/donor/site), whether the downstream use needs a probability or a label, and whether you are building a prediction model (this skill) or testing a treatment effect (clinical-biostatistics).

## Quick Start

Tell your AI agent what you want to do:
- "Run nested cross-validation and report AUC with an interval"
- "Keep all samples from each patient in the same fold"
- "Check whether my predicted probabilities are calibrated, not just discriminating"
- "Is my high-AUC model actually clinically useful? Run a decision curve"

## Example Prompts

### Nested CV and leakage

> "Run nested cross-validation on my expression classifier with feature selection and scaling inside the folds, 5 outer folds and 3 inner folds, and report AUC mean and spread."

> "Audit my pipeline for leakage: am I fitting any normalization, batch correction, or feature selection before the split?"

### Matching the question

> "My dataset has multiple biopsies per patient. Use a group-aware split so no patient spans folds."

> "I want to know if my model transports to a new hospital. Set up leave-one-site-out validation."

### Calibration and utility

> "Plot a reliability curve and compute the Brier score for my classifier, then recalibrate on a held-out fold."

> "My model has AUC 0.85 but I need risk estimates. Check calibration and run decision-curve analysis for net benefit."

### Imbalance and reporting

> "My positive class is rare. Should I use SMOTE, and which metrics should I report?"

> "Prepare the validation reporting to TRIPOD+AI: discrimination, calibration, subgroup performance, and uncertainty."

## What the Agent Will Do

1. Choose the CV scheme from the generalization question (random, group, site, temporal)
2. Wrap all preprocessing and feature selection in a `Pipeline` so they refit per fold
3. Use nested CV whenever hyperparameters or options are tuned
4. Report discrimination with an interval, plus calibration (reliability curve, Brier) and, for decisions, net benefit
5. Flag leakage, threshold-on-test, and resampling-for-imbalance pitfalls
6. Point to external/temporal validation and TRIPOD+AI reporting for clinical claims

## Tips

- A clean train/test split does not prevent leakage if a scaler, ComBat, or a duplicate patient already contaminated the test set
- AUC is invariant to monotone transforms of the score, so it cannot tell you whether the probabilities are honest -- check calibration separately
- The calibration slope is the single most informative calibration number; <1 means overfitting
- SMOTE/oversampling destroys calibration for no AUC gain on risk models; move the threshold on a calibrated model instead
- Choosing the classification threshold on the test set is leakage; pick it on a separate fold and report the locked threshold once
- "10 events per variable" is obsolete; size the sample with the Riley framework (`pmsampsize`)
- TRIPOD+AI (2024) is the reporting target; it requires calibration and fairness, not just discrimination

## Related Skills

- machine-learning/biomarker-discovery - Feature selection run inside the CV fold
- machine-learning/omics-classifiers - Model training, calibration directions, and imbalance handling
- machine-learning/survival-analysis - Validation metrics for time-to-event models
- experimental-design/batch-design - Designing out batch-outcome confounding before analysis
- experimental-design/multiple-testing - FDR control for high-dimensional testing
- clinical-biostatistics/trial-reporting - Confirmatory-trial reporting and the prediction-vs-inference boundary

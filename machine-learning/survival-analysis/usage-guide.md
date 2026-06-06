# Predictive Survival Modeling Usage Guide

## Overview

Build and validate individualized time-to-event risk models on clinical and omics data with penalized Cox, random survival forests, gradient-boosted and deep survival models. The discipline this skill enforces is that the C-index alone is the cardinal evaluation sin: it is censoring-dependent, blind to calibration, and insensitive, so a credible model is judged by Uno's C plus time-dependent AUC, integrated Brier versus a Kaplan-Meier baseline, and calibration. For Kaplan-Meier, log-rank, and classical Cox hazard-ratio inference, use clinical-biostatistics/survival-analysis.

## Prerequisites

```bash
pip install scikit-survival lifelines pandas
# pip install pycox torch   # only for deep survival models (DeepSurv, DeepHit)
```

Conceptual prerequisites: a time and a boolean event column; awareness of competing events (which require the cumulative incidence function, not 1-KM) and of immortal time bias (no group defined by a post-baseline event); and feature selection placed inside the resampling loop.

## Quick Start

Tell your AI agent what you want to do:
- "Fit a penalized Cox prognostic signature and validate it properly"
- "Compare a random survival forest against an elastic-net Cox baseline"
- "Evaluate my survival model with Uno's C, time-dependent AUC, and integrated Brier, not just the C-index"
- "I have competing risks -- model the cumulative incidence correctly"

## Example Prompts

### Building predictors

> "Build an elastic-net Cox prognostic signature from my expression matrix and clinical covariates, with feature selection inside nested cross-validation."

> "Fit a random survival forest and a penalized Cox model and tell me whether the forest actually beats the linear baseline on held-out data."

### Prediction-grade evaluation

> "Evaluate my survival model with Uno's IPCW C at a 5-year horizon, time-dependent AUC over the follow-up, integrated Brier versus a Kaplan-Meier baseline, and a calibration curve."

> "My model has a high C-index. Check whether its predicted absolute risks are actually calibrated."

### Competing risks and bias

> "Patients also die of other causes. Model the cumulative incidence with Fine-Gray and report cause-specific hazards too, and do not use 1 minus Kaplan-Meier."

> "My exposure is defined by something that happens after baseline. Check for immortal time bias and fix it with landmarking."

## What the Agent Will Do

1. Build the structured target (boolean event, float time) and choose a model from the taxonomy
2. Fit penalized Cox and RSF baselines before escalating to deep models
3. Keep feature selection and tuning inside nested CV with survival metrics
4. Evaluate with Uno's C(tau), time-dependent AUC, IBS versus KM, and calibration at clinical horizons
5. Handle competing risks with the CIF and report both cause-specific and Fine-Gray models
6. Check censoring assumptions, immortal time bias, and use landmarking for dynamic prediction

## Tips

- The C-index is invariant to any monotone transform of the risk score, so it cannot tell you whether the absolute risks are right -- always check calibration
- Harrell's C is biased upward under heavy censoring; report Uno's IPCW C with an explicit truncation time
- On typical clinical-omics n, penalized Cox and RSF are very hard to beat; benchmark deep models against them
- Kaplan-Meier overestimates incidence under competing risks; report the cumulative incidence function
- A Fine-Gray subdistribution hazard ratio is not the cause-specific hazard ratio and not an effect on the event rate
- Internal time-varying covariates cannot be plugged into the future; use landmarking for dynamic prediction
- scikit-survival needs a structured array (boolean event, float time) and the training y first for IPCW metrics

## Related Skills

- clinical-biostatistics/survival-analysis - Kaplan-Meier, log-rank, classical Cox inference, PH diagnostics for trials
- machine-learning/model-validation - Nested CV, calibration, and optimism correction shared with prediction models
- machine-learning/biomarker-discovery - Selecting prognostic genes inside the resampling loop
- differential-expression/de-results - Pre-filter candidate prognostic genes
- clinical-databases/variant-prioritization - Clinical interpretation of prognostic variants

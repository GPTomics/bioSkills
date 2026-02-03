# Survival Prediction Usage Guide

## Overview

Analyze time-to-event data using Kaplan-Meier survival curves, log-rank tests for group comparisons, and Cox proportional hazards regression for multivariate survival modeling.

## Prerequisites

```bash
pip install lifelines pandas matplotlib
```

## Quick Start

Tell your AI agent what you want to do:
- "Plot Kaplan-Meier survival curves for my patient groups"
- "Compare survival between high and low risk groups with log-rank test"
- "Build a Cox regression model from clinical and expression data"
- "Calculate the concordance index for my survival model"

## Example Prompts

### Kaplan-Meier Analysis

> "I have clinical data with survival_time and event columns. Plot Kaplan-Meier curves comparing treatment vs control groups and calculate the log-rank p-value."

> "Create survival curves for my high and low risk groups. Add median survival times to the plot."

### Cox Regression

> "Build a Cox proportional hazards model using age, stage, and the expression of BRCA1 and TP53. Show me the hazard ratios."

> "Run a multivariate Cox model on my clinical data. Which features significantly predict survival?"

### Prognostic Signatures

> "Screen all genes in my expression matrix for survival association. Which genes have hazard ratio > 2 and p < 0.01?"

### Risk Stratification

> "Calculate Cox risk scores for my patients and split into high/low risk groups at the median. Show the KM plot."

## What the Agent Will Do

1. Load clinical data with survival time and event status
2. Fit Kaplan-Meier estimator per group
3. Perform log-rank test for group comparisons
4. Build Cox regression model with covariates
5. Report hazard ratios, confidence intervals, p-values
6. Evaluate model with concordance index

## Tips

- Event=1 means the event occurred; Event=0 means censored (no event by end of follow-up)
- Log-rank test is non-parametric; Cox is semi-parametric with proportional hazards assumption
- Check proportional hazards assumption with `cph.check_assumptions()`
- Use penalizer (L2 regularization) when features outnumber samples
- C-index > 0.7 is generally considered useful for prediction
- Note: lifelines `sklearn_adaptor` was removed in v0.28; use native API

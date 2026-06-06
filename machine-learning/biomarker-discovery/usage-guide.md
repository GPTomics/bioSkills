# Feature Selection for Biomarker Discovery Usage Guide

## Overview

Select candidate biomarker features from high-dimensional omics data using Boruta (all-relevant), mRMR, LASSO/elastic-net (minimal-optimal), and stability selection. The discipline this skill enforces is that most gene signatures do not replicate: selection must happen inside the cross-validation loop, the right null is a random size-matched signature rather than "no association," and a selected list is reported with a stability index, never on its own.

## Prerequisites

```bash
pip install Boruta mrmr-selection scikit-learn pandas numpy
```

Conceptual prerequisites: decide whether you want all-relevant (which genes carry signal) or minimal-optimal (smallest predictive set) before choosing a method; for single-cell data, pseudobulk to the donor level first.

## Quick Start

Tell your AI agent what you want to do:
- "Select all-relevant biomarker genes with Boruta for pathway interpretation"
- "Build a small stable elastic-net signature and report its stability index"
- "Run leakage-safe selection inside cross-validation and give an honest AUC"
- "Is my selected gene list reproducible, or a resampling accident?"

## Example Prompts

### All-relevant vs minimal-optimal

> "I want to enumerate every gene implicated in disease status for downstream pathway analysis, including redundant co-expressed genes. Which selector should I use?"

> "I need the smallest possible panel for an assay. Use elastic net rather than bare LASSO and tell me why."

### Leakage-safe estimation

> "Select the top genes and estimate classifier performance without leakage. Put selection inside the cross-validation folds and use nested CV for tuning."

### Reproducibility

> "Run subsampling stability selection, keep features selected in more than 60% of runs, and report the Nogueira stability index next to the accuracy."

> "My signature barely overlaps a published one. Is that a problem?"

### Correct nulls

> "Does my signature beat a size-matched random signature and a proliferation meta-gene in independent data?"

## What the Agent Will Do

1. Establish the question: all-relevant (Boruta) vs minimal-optimal (elastic net) vs stable signature
2. For single-cell, pseudobulk to the donor level
3. Run selection inside a `Pipeline` so held-out folds never inform feature choice
4. Estimate performance by nested CV; report discrimination with an interval
5. Quantify stability across subsamples (selection frequency + a chance-corrected index)
6. Benchmark against a random-signature null and recommend independent-cohort validation

## Tips

- Absence from a LASSO/minimal-optimal set is not evidence a gene is irrelevant -- it may be redundant with a kept gene
- "LASSO kept gene X but not its co-expressed partner" is L1 geometry, not biology; use elastic net or report selection frequencies
- Selecting features on the whole dataset before CV inflates AUC to near-perfect even on pure noise -- selection must be inside the fold
- A signature being "significantly associated with outcome" is weak evidence; random signatures clear that bar (Venet 2011)
- Report a stability index (Nogueira 2018) alongside accuracy; a high-accuracy, low-stability signature is an accident
- Effect sizes at selected features are inflated (winner's curse); estimate them on an independent split and size replication for the shrunken effect
- BorutaPy needs numpy arrays and a tree estimator; it returns a redundant all-relevant set by design

## Related Skills

- machine-learning/model-validation - Nested CV and leakage-safe estimation of the selected model
- machine-learning/prediction-explanation - Why SHAP rankings are not a validated selection method
- machine-learning/omics-classifiers - Build a classifier from the selected features
- differential-expression/de-results - Pre-filter candidates with differential expression
- experimental-design/multiple-testing - FDR control and why it is orthogonal to selection stability
- experimental-design/power-analysis - Sample size for a stable signature vs an accurate predictor
- pathway-analysis/go-enrichment - Functional enrichment of an all-relevant gene set

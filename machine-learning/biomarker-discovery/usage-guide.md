# Feature Selection Usage Guide

## Overview

Select informative features for biomarker discovery using Boruta all-relevant selection, mRMR minimum redundancy, and LASSO regularization.

## Prerequisites

```bash
pip install Boruta mrmr-selection scikit-learn pandas numpy
```

## Quick Start

Tell your AI agent what you want to do:
- "Select biomarker genes using Boruta"
- "Find a minimal set of non-redundant features with mRMR"
- "Use LASSO to find sparse biomarkers"
- "Which features are stable across bootstrap samples?"

## Example Prompts

### Boruta Selection

> "Run Boruta on my expression matrix to find all genes relevant for predicting disease status. I want the complete set of biomarkers, not just a minimal subset."

> "Use Boruta feature selection on my data. How many genes are selected vs tentative?"

### mRMR Selection

> "Select 50 non-redundant biomarker genes using mRMR from my expression data."

### LASSO Biomarkers

> "Use LASSO with cross-validation to find a sparse set of genes that predict my outcome. What alpha was selected?"

### Stability Selection

> "Run stability selection with 100 bootstrap samples. Which features are selected in more than 60% of runs?"

### Combined Approaches

> "First filter to top 5000 genes by ANOVA, then run Boruta on the filtered set."

## What the Agent Will Do

1. Load expression matrix and sample labels
2. Apply appropriate feature selection method
3. Report number of selected features
4. Rank features by importance/stability
5. Save selected feature list

## Tips

- Boruta finds ALL relevant features; mRMR/LASSO find minimal sets
- Pre-filter with univariate tests (top 1000-5000) before Boruta on large matrices
- mRMR is good when you want exactly K features with low redundancy
- LASSO picks arbitrarily among correlated features; stability selection helps
- Consider running multiple methods and taking intersection or union
- Validate selected features on independent data or nested CV

# Statistical Analysis - Usage Guide

## Overview

Statistical analysis identifies metabolites associated with biological conditions. Methods range from simple univariate tests to multivariate models for biomarker discovery.

## Prerequisites

```bash
# Python
pip install scipy statsmodels numpy pandas matplotlib scikit-learn

# R packages
install.packages(c("mixOmics", "ropls", "pROC"))
```

## Quick Start

Tell your AI agent what you want to do:
- "Find differentially abundant metabolites between treatment and control"
- "Run PLS-DA and identify important features by VIP score"

## Example Prompts

### Univariate Analysis
> "Run Welch's t-tests on my log2-transformed metabolomics data with BH FDR correction"
> "Perform ANOVA across my three treatment groups and identify significant metabolites"
> "Calculate log2 fold changes and create a volcano plot"

### Multivariate Analysis
> "Run PCA for exploratory analysis and check sample grouping"
> "Build a PLS-DA model with 10-fold cross-validation and calculate VIP scores"
> "Use OPLS-DA for biomarker discovery between disease and healthy groups"

### Biomarker Selection
> "Identify metabolites with VIP > 1, FDR < 0.05, and |log2FC| > 1"
> "Calculate ROC curves and AUC for top candidate biomarkers"
> "Build a Random Forest classifier and rank feature importance"

### Model Validation
> "Validate my PLS-DA model with 100 permutation tests"
> "Report Q2 and R2 from cross-validation"

## What the Agent Will Do

1. Check data distribution and apply appropriate tests
2. Run univariate tests with multiple testing correction
3. Build and validate multivariate models
4. Calculate VIP scores and feature importance
5. Generate summary tables and visualizations
6. Export significant metabolites

## Tips

- Log2-transform raw intensities before statistical testing; compute fold change as difference of means on the log2 scale (not log2 of ratio of raw means)
- Always use Welch's t-test over Student's; in Python, set `equal_var=False` in `scipy.stats.ttest_ind()`
- Always correct for multiple testing (FDR/BH method is standard)
- Validate PLS-DA with permutation testing (Q2 should exceed permuted values)
- VIP > 1 is common threshold, but combine with FDR for confidence
- Use 5-10 fold CV with 50+ repeats for stable model assessment
- Report both univariate (FDR) and multivariate (VIP) evidence

## Method Selection

| Method | Samples | Use Case |
|--------|---------|----------|
| Welch's t-test | 2 groups | Default for two-group comparison |
| limma moderated t | 2 groups, small n | Better with < 5 samples per group |
| Mann-Whitney U | 2 groups | Non-parametric alternative |
| Welch's ANOVA | 3+ groups | Multiple conditions |
| PCA | Any | Exploratory, QC |
| PLS-DA | 2+ groups | Classification, VIP |
| OPLS-DA | 2 groups | Biomarker discovery |
| Random Forest | 2+ groups | Non-linear importance |

## Biomarker Criteria

| Metric | Threshold |
|--------|-----------|
| FDR | < 0.05 (discovery: < 0.1) |
| log2FC | > 1 (2-fold change) |
| VIP | > 1 |
| AUC | > 0.7 (moderate), > 0.8 (good) |

## References

- mixOmics: doi:10.1371/journal.pcbi.1005752
- ropls (OPLS-DA): doi:10.1021/acs.jproteome.5b00354

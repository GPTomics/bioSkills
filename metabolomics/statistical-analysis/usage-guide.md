# Statistical Analysis Usage Guide

## Overview

Statistical analysis identifies metabolites associated with biological conditions. Methods range from simple univariate tests to complex multivariate models.

## Method Selection

| Method | Samples | Use Case |
|--------|---------|----------|
| t-test | 2 groups | Simple comparison |
| ANOVA | 3+ groups | Multiple conditions |
| PCA | Any | Exploratory, QC |
| PLS-DA | 2+ groups | Classification, VIP |
| OPLS-DA | 2 groups | Biomarker discovery |
| Random Forest | 2+ groups | Non-linear, importance |

## Multiple Testing

Always correct for multiple testing:
- **FDR (BH)**: Most common, controls false discovery rate
- **Bonferroni**: Conservative, controls family-wise error
- **q-value**: Similar to FDR

## VIP Scores

Variable Importance in Projection (VIP):
- VIP > 1: Important (common threshold)
- VIP > 1.5: Very important
- Use with FDR for confidence

## Cross-Validation

Always validate models:
- 5-10 fold CV
- 50+ repeats for stability
- Report Q2 (predictive ability)

## Biomarker Criteria

Typical thresholds:
- FDR < 0.05 (or 0.1 for discovery)
- |log2FC| > 1 (2-fold change)
- VIP > 1
- AUC > 0.7 (moderate), > 0.8 (good)

## References

- mixOmics: doi:10.1371/journal.pcbi.1005752
- ropls (OPLS-DA): doi:10.1021/acs.jproteome.5b00354

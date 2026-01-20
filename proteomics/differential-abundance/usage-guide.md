# Differential Protein Abundance

## Overview

Statistical testing identifies proteins with significantly different abundance between experimental conditions.

## Workflow

1. **Normalized, imputed intensities** from quantification
2. **Define contrasts** (which groups to compare)
3. **Statistical testing** (t-test, limma, linear models)
4. **Multiple testing correction** (Benjamini-Hochberg FDR)
5. **Apply thresholds** (adj. p-value < 0.05, |log2FC| > 1)

## Statistical Methods

| Method | Description | Best For |
|--------|-------------|----------|
| t-test | Simple two-group comparison | Few replicates |
| limma | Empirical Bayes moderated t-test | Small sample sizes |
| MSstats | Mixed-effects models | Complex designs |
| DEP | Full proteomics workflow | Beginners |

## Multiple Testing Correction

- **Bonferroni**: Conservative, controls FWER
- **Benjamini-Hochberg**: Controls FDR, more powerful
- **q-value**: Storey's method, permutation-based

## Significance Thresholds

Typical thresholds for proteomics:
- **Adjusted p-value**: < 0.05 (or 0.01 for stringent)
- **Log2 fold change**: > 1 (2-fold change) or > 0.58 (1.5-fold)

## Common Pitfalls

1. **Not accounting for missing values**: Can bias results
2. **Ignoring batch effects**: Include batch in model
3. **Too few replicates**: Need n >= 3 per group for statistics
4. **Multiple testing**: Always correct p-values
5. **Outliers**: Check with PCA, consider robust methods

## Reporting

Include:
- Number of proteins tested
- Significance thresholds used
- Number of up/down regulated
- Normalization and imputation methods
- Software versions

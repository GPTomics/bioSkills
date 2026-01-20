# Differential Abundance - Usage Guide

## Overview
Identify proteins with significantly different abundance between experimental conditions using statistical testing and multiple testing correction.

## Prerequisites
```bash
pip install numpy pandas scipy statsmodels
# R packages: BiocManager::install(c("limma", "MSstats", "DEP", "proDA", "QFeatures"))
```

## Quick Start
Tell your AI agent what you want to do:
- "Find differentially abundant proteins between treatment and control"
- "Run limma analysis on my normalized protein matrix"
- "Generate a volcano plot showing significant proteins"

## Example Prompts

### Statistical Testing
> "Run limma differential analysis comparing treatment vs control groups"

> "Perform a t-test on each protein with Benjamini-Hochberg correction"

> "Use MSstats to test for differential abundance in my label-free data"

### Complex Designs
> "Set up a limma model with treatment and batch as covariates"

> "Test for interaction effects between drug and timepoint"

> "Run paired differential analysis for my before/after samples"

### Results Filtering
> "Filter to proteins with adjusted p-value < 0.05 and |log2FC| > 1"

> "Extract the top 50 most significantly changed proteins"

> "Identify proteins changing in the same direction across all comparisons"

### Visualization
> "Create a volcano plot highlighting the significant proteins"

> "Make a heatmap of the top differentially abundant proteins"

## What the Agent Will Do
1. Load normalized, imputed protein matrix
2. Define experimental design and contrasts
3. Fit statistical model (limma/t-test/MSstats)
4. Apply multiple testing correction (BH FDR)
5. Filter by significance thresholds
6. Generate results table and visualizations

## Statistical Methods

| Method | Description | Best For |
|--------|-------------|----------|
| t-test | Simple two-group comparison | Large sample sizes |
| limma | Empirical Bayes moderated t-test | Small sample sizes |
| MSstats | Mixed-effects models | Complex designs |
| proDA | Probabilistic dropout model | High missing values |

## Significance Thresholds

Typical thresholds for proteomics:
- **Adjusted p-value**: < 0.05 (or 0.01 for stringent)
- **Log2 fold change**: > 1 (2-fold) or > 0.58 (1.5-fold)

## Tips
- Always use adjusted p-values (not raw p-values)
- Include batch as covariate if samples were processed separately
- Need n >= 3 replicates per group for reliable statistics
- Check volcano plot for bias (symmetric up/down regulation expected)
- Report number tested, thresholds, and number significant

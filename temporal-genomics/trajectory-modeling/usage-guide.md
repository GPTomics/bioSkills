# Temporal Trajectory Modeling - Usage Guide

## Overview

Fits continuous non-linear curves to gene expression time series using generalized additive models (GAMs) and detects abrupt shifts with changepoint algorithms. Enables formal statistical testing of temporal trends, condition-dependent trajectory differences, and dynamic breakpoints in temporal expression data.

## Prerequisites

### R
```r
install.packages(c('mgcv', 'segmented'))
BiocManager::install('tradeSeq')
```

### Python
```bash
pip install ruptures numpy matplotlib
```

### Data Requirements
- Expression values with corresponding timepoint metadata
- For GAMs: at least 6 unique timepoints (more is better for smooth estimation)
- For condition comparison: expression data with both timepoint and condition labels
- For changepoint detection: ordered time-series expression values

## Quick Start

Tell the AI agent what to model:
- "Fit smooth curves to gene expression over time using GAMs"
- "Test whether my treated samples have a different temporal trajectory than controls"
- "Detect changepoints in temporal gene expression where dynamics shift abruptly"
- "Compare GAM vs linear model fits for my time-course genes"

## Example Prompts

### GAM Trajectory Fitting
> "Fit a generalized additive model to each gene in my RNA-seq time-course data and test which genes have significant temporal trends."

> "I have 10 timepoints of expression data. Fit smooth curves and give me the effective degrees of freedom for each gene."

### Condition Comparison
> "I have time-course RNA-seq data from WT and KO mice. Test which genes have different temporal trajectories between genotypes."

> "Compare the temporal expression curves of drug-treated vs untreated samples using GAMs with difference smooths."

### Changepoint Detection
> "Find genes where expression dynamics change abruptly during my developmental time course."

> "Detect breakpoints in my temporal gene expression data -- when does the expression regime shift?"

### Model Selection
> "For each gene, compare linear, quadratic, and GAM fits to my time-series data and report which model fits best by AIC."

## What the Agent Will Do

1. Load expression data and timepoint metadata
2. Fit GAM or spline models per gene with appropriate basis dimensions
3. Run model diagnostics (gam.check, concurvity, residual plots)
4. Test for significant temporal trends or condition differences
5. Correct p-values for multiple testing (BH FDR)
6. Detect changepoints if abrupt shifts are expected
7. Generate trajectory plots with confidence intervals
8. Export results tables with test statistics and model parameters

## Tips

- Use REML estimation for GAMs (method='REML') rather than GCV; REML is more robust to overfitting
- The basis dimension k should be less than the number of unique timepoints; k=6 works well for 8-12 timepoints
- Always check gam.check() output: if k-index < 1.0, the smooth may be too rigid and k should be increased
- For condition comparison, the difference smooth approach is more interpretable than fitting fully separate models
- Changepoint detection is most useful for developmental or treatment time courses with expected regime shifts
- Combine GAMs (smooth trends) with changepoint detection (abrupt shifts) for comprehensive trajectory analysis
- BIC-based penalty in ruptures Pelt is standard; increase penalty for fewer, more confident changepoints

## Related Skills

- temporal-genomics/temporal-clustering - Group genes after trajectory fitting
- temporal-genomics/circadian-rhythms - Periodic trajectory models
- differential-expression/timeseries-de - Linear model alternatives for temporal DE
- single-cell/trajectory-inference - Single-cell pseudotime trajectories

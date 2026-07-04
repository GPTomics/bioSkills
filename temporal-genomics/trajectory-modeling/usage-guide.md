# Temporal Trajectory Modeling - Usage Guide

## Overview

This skill fits continuous trajectories to BULK or time-resolved omics where the x-axis is MEASURED experimental time (hours, days, developmental stage) or a pseudobulk value aggregated over real time - it is explicitly not single-cell pseudotime. Two complementary jobs: penalized-spline GAMs (mgcv) model gradually curving trends and test whether expression changes or differs between conditions, while changepoint methods (segmented, ruptures) locate discrete regime shifts. Three decisions gate correctness before any p-value: the distributional scale (raw overdispersed counts need `nb()` plus a library-size offset, not Gaussian), residual autocorrelation across ordered timepoints (which inflates smooth-term significance if ignored), and whether the mechanism is gradual (GAM) or genuinely abrupt (changepoint). The basis dimension k is a flexibility ceiling below the number of timepoints, not a knot count; REML picks the realized wiggliness and edf reports it.

## Prerequisites

### R
```r
install.packages(c('mgcv', 'segmented'))
BiocManager::install('tradeSeq')
```

### Python
```bash
pip install ruptures numpy matplotlib pandas
```

### Data Requirements
- Expression values with corresponding measured-time metadata (real time, not inferred pseudotime)
- For GAMs: at least 6 unique timepoints (more is better); k must stay below the number of unique timepoints
- Raw counts require the negative-binomial family with a library-size offset; a Gaussian GAM needs variance-stabilized values (vst, rlog, log-CPM)
- For condition comparison: expression data with both timepoint and condition labels
- For changepoint detection: ordered time-series expression values

## Quick Start

Tell the AI agent what to model:
- "Fit smooth GAM curves to my gene expression over time and test which genes change"
- "My input is raw counts - model them on the correct distributional scale"
- "Test whether treated samples diverge from controls over time with a difference smooth"
- "My timepoints are correlated - handle residual autocorrelation in the GAM"
- "Detect changepoints where expression shifts abruptly, and check the shift is real"

## Example Prompts

### GAM Trajectory Fitting
> "Fit a generalized additive model to each gene in my RNA-seq time-course and report which genes have significant temporal trends after FDR correction."

> "I have 10 timepoints of variance-stabilized expression. Fit smooth curves and give me the effective degrees of freedom per gene, using REML."

### Distribution and Autocorrelation
> "My data are raw RNA-seq counts with varying library sizes. Fit a negative-binomial GAM with a library-size offset rather than a Gaussian model."

> "My samples are a single series measured repeatedly over time, so timepoints are correlated. Fit a GAM that models the residual autocorrelation and check the residual ACF."

### Condition Comparison
> "I have time-course RNA-seq from WT and KO mice. Test which genes have diverging temporal trajectories using an ordered-factor difference smooth."

> "Compare drug-treated vs untreated expression curves over time and give me a single p-value for whether the trajectories differ."

### Changepoint Detection
> "Find genes where expression shifts abruptly during my developmental time course, and confirm the shift is a real regime change rather than a smooth ramp."

> "Detect breakpoints in my temporal expression data. Distinguish a slope change from a level jump, and tell me when the regime shifts."

### Model Selection
> "For each gene, compare a smooth GAM against a linear fit by AIC and report whether non-linear modeling is warranted."

## What the Agent Will Do

1. Load expression data and measured-time metadata; confirm the scale matches the family (vst/log-CPM for Gaussian, raw counts with an offset for nb())
2. Fit GAMs per gene with k set below the number of unique timepoints and REML smoothing
3. Model residual autocorrelation with corAR1/bam(rho=) where the design warrants it, and inspect the residual ACF
4. Run diagnostics (gam.check/k.check, concurvity for multi-smooth models), responding to a low k-index by doubling k and refitting rather than reflexively raising it
5. Test condition divergence with an ordered-factor difference smooth plus the required parametric main effect
6. Correct smooth-term p-values across genes with BH FDR, treating them as approximate
7. Detect changepoints with a noise-calibrated penalty and a cost model matched to the shift type, pre-testing that a break exists
8. Predict trajectories only within the sampled time range with pointwise confidence bands, and export result tables

## Tips

- Set k as a ceiling (generously, but below the number of unique timepoints) and let REML choose the wiggliness; read edf, not k, for realized complexity
- Never fit a Gaussian GAM on raw counts - use `family=nb()` with `offset(log(library_size))`, or model vst/log-CPM values
- Ignoring residual autocorrelation inflates smooth-term p-values; model it or at least check the lag-1 residual ACF
- A low k-index does not automatically mean k is too low - double k and refit; if edf barely moves, suspect autocorrelation or the distribution
- Smooth-term p-values are approximate; treat them as categorical significant/not and apply FDR across genes
- For condition comparison, an ordered-factor difference smooth gives a single divergence p-value; keep the parametric main effect
- In ruptures, estimate noise from `np.var(np.diff(signal))/2`, not total variance, or the BIC penalty under-detects; use `model='l2'` for level shifts (BIC is derived for it) and calibrate `'rbf'` empirically
- Pre-test that a changepoint exists (davies.test) or compare piecewise vs smooth-GAM AIC before interpreting a break - a liberal penalty always finds breaks in smooth data
- Never predict or interpret outside `[min(time), max(time)]`; the spline and its interval diverge past the data

## Related Skills

- temporal-genomics/temporal-clustering - group genes by trajectory shape after fitting
- temporal-genomics/circadian-rhythms - periodic (known-period) trajectory models rather than smooth trends
- temporal-genomics/periodicity-detection - discover unknown-period oscillation instead of a smooth trend
- differential-expression/timeseries-de - linear/spline model alternatives for temporal DE
- single-cell/trajectory-inference - single-cell pseudotime, the case tradeSeq is built for

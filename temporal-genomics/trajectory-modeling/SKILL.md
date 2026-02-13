---
name: bio-temporal-genomics-trajectory-modeling
description: Models continuous temporal trajectories from bulk or time-resolved omics data using generalized additive models (mgcv), spline regression, and changepoint detection (segmented, ruptures). Fits smooth gene expression curves and tests trajectory differences between conditions. Use when fitting non-linear temporal models to bulk time-series data or comparing developmental trajectories across conditions. Not for single-cell pseudotime (see single-cell/trajectory-inference).
tool_type: mixed
primary_tool: mgcv
---

## Version Compatibility

Reference examples tested with: R stats (base), numpy 1.26+, pandas 2.2+, scanpy 1.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Temporal Trajectory Modeling

**"Fit smooth curves to my gene expression time series"** → Model continuous temporal trajectories using generalized additive models (GAMs) or spline regression, test for condition differences, and detect changepoints where dynamics shift abruptly.
- R: `mgcv::gam()` for GAM fitting with smooth terms
- Python: `ruptures` for changepoint detection in temporal profiles

Fits smooth non-linear curves to gene expression time series using generalized additive models (GAMs) and detects abrupt changes in temporal dynamics using changepoint algorithms.

## Core Workflow

1. Prepare expression data with timepoint and condition metadata
2. Fit GAM or spline models per gene
3. Test for significant temporal trends and condition differences
4. Detect changepoints where trajectory dynamics shift
5. Predict and visualize fitted trajectories with confidence intervals

## mgcv GAM (R)

**Goal:** Fit smooth non-linear curves to gene expression time series and test for significant temporal trends or condition differences.

**Approach:** Use generalized additive models with penalized smooth terms to capture non-linear dynamics, compare trajectories between conditions using interaction smooths, and extract predicted values with confidence intervals.

### Basic GAM Fitting (R stats (base)+)

```r
library(mgcv)

# gam() with s() smooth terms for non-linear temporal dynamics
# k=6: basis dimension (number of knots); controls smoothness
# k should be less than the number of unique timepoints
# Too low k = underfit; too high k = overfit; k=6 is a good default for 8-12 timepoints
fit <- gam(expression ~ s(time, k = 6), data = gene_df, method = 'REML')

summary(fit)
# edf (effective degrees of freedom): edf near 1 = linear; edf near k-1 = highly non-linear
# p-value for s(time): tests whether the smooth term is significantly non-zero
```

### Condition Comparison with GAM (R stats (base)+)

```r
# 'by' argument fits separate smooths per condition
# This tests whether temporal trajectories differ between groups
fit_cond <- gam(
    expression ~ condition + s(time, k = 6, by = condition),
    data = gene_df, method = 'REML'
)

# Difference smooth: directly models the trajectory difference
fit_diff <- gam(
    expression ~ s(time, k = 6) + s(time, k = 6, by = is_treated),
    data = gene_df, method = 'REML'
)

# The second smooth (by=is_treated) represents the deviation from baseline
# Its p-value tests whether conditions differ in their temporal trajectory
summary(fit_diff)
```

### Model Diagnostics (R stats (base)+)

```r
# gam.check() provides residual plots and basis dimension checks
# k-index: if < 1.0, consider increasing k (smooth may be too rigid)
gam.check(fit)

# Concurvity: analogous to collinearity for smooth terms
# Values > 0.8 indicate redundancy between smooth terms
concurvity(fit_cond, full = TRUE)
```

### Prediction and Visualization (R stats (base)+)

```r
# Predict on fine time grid with confidence intervals
new_data <- data.frame(time = seq(min(gene_df$time), max(gene_df$time), length.out = 200))
pred <- predict(fit, newdata = new_data, se.fit = TRUE)

new_data$fitted <- pred$fit
# 1.96 * SE: 95% confidence interval assuming approximate normality
new_data$lower <- pred$fit - 1.96 * pred$se.fit
new_data$upper <- pred$fit + 1.96 * pred$se.fit
```

### Genome-Wide GAM Fitting (R stats (base)+)

```r
# Fit GAMs for all genes and test for significant temporal trends
results <- data.frame()
for (gene in rownames(expr_mat)) {
    gene_df <- data.frame(expression = as.numeric(expr_mat[gene, ]), time = timepoints)
    fit <- gam(expression ~ s(time, k = 6), data = gene_df, method = 'REML')
    s_table <- summary(fit)$s.table
    results <- rbind(results, data.frame(
        gene = gene, edf = s_table[, 'edf'],
        F_stat = s_table[, 'F'], p_value = s_table[, 'p-value']
    ))
}

results$q_value <- p.adjust(results$p_value, method = 'BH')
# q < 0.05: standard FDR threshold for temporal significance
temporal_genes <- results[results$q_value < 0.05, ]
```

## tradeSeq (R/Bioconductor)

Wrapper around mgcv designed for trajectory analysis with built-in statistical tests.

```r
library(tradeSeq)

# fitGAM for pseudobulk or aggregated time-course data
# nKnots=6: number of interior knots; similar role to k in mgcv
# nKnots should be less than the number of unique timepoints
sce <- fitGAM(counts = count_mat, pseudotime = time_mat, cellWeights = weight_mat, nKnots = 6)

# associationTest: tests whether expression changes over time (any gene)
assoc_res <- associationTest(sce)

# conditionTest: tests trajectory differences between conditions
cond_res <- conditionTest(sce)
```

## segmented (R)

Piecewise linear regression with automatic breakpoint detection.

```r
library(segmented)

# Fit linear model first, then test for breakpoints
lm_fit <- lm(expression ~ time, data = gene_df)

# segmented() searches for a breakpoint in the linear relationship
# psi: initial guess for breakpoint location; NA lets it estimate
seg_fit <- segmented(lm_fit, seg.Z = ~time, psi = NA)

# davies.test: tests whether a breakpoint exists (H0: no breakpoint)
# p < 0.05: evidence of a changepoint in the temporal trajectory
davies.test(lm_fit, seg.Z = ~time)

# Extract breakpoint estimate and confidence interval
summary(seg_fit)$psi
```

## ruptures (Python)

Changepoint detection for identifying abrupt shifts in temporal dynamics.

```python
import numpy as np
import ruptures as rpt

# signal: 1D array of gene expression over time
signal = np.array(expression_values)

# Pelt: penalized exact linear time algorithm
# model='rbf': radial basis function kernel; detects mean and variance changes
# min_size=2: minimum segment length; at least 2 timepoints per segment
algo = rpt.Pelt(model='rbf', min_size=2).fit(signal)

# pen=np.log(n) * sigma**2: BIC penalty; balances fit vs complexity
# Higher penalty = fewer changepoints (more conservative)
# BIC-based penalty is standard for model selection in changepoint detection
n = len(signal)
penalty = np.log(n) * np.var(signal)
changepoints = algo.predict(pen=penalty)
```

### Binary Segmentation Alternative (R stats (base)+)

```python
# BinSeg: faster approximate method for long series
# n_bkps=3: maximum number of changepoints to detect
# Use when computational speed matters or as a sanity check against Pelt
algo_binseg = rpt.Binseg(model='rbf', min_size=2).fit(signal)
changepoints_binseg = algo_binseg.predict(n_bkps=3)
```

### Genome-Wide Changepoint Detection

```python
import pandas as pd

results = []
for gene_idx in range(expr_mat.shape[0]):
    signal = expr_mat[gene_idx, :]
    algo = rpt.Pelt(model='rbf', min_size=2).fit(signal)
    penalty = np.log(len(signal)) * np.var(signal)
    bkps = algo.predict(pen=penalty)
    # bkps[-1] is always n (end of signal); preceding values are changepoint indices
    n_changes = len(bkps) - 1
    results.append({'gene': f'gene_{gene_idx}', 'n_changepoints': n_changes,
                     'changepoint_indices': bkps[:-1]})

results_df = pd.DataFrame(results)
```

## Model Comparison

| Method | Model Type | Best For | Key Parameter |
|--------|-----------|----------|---------------|
| mgcv GAM | Smooth non-linear | Continuous trajectories | k (basis dimension) |
| tradeSeq | GAM wrapper | Condition comparison | nKnots |
| segmented | Piecewise linear | Breakpoint detection | psi (initial guess) |
| ruptures Pelt | Changepoint | Abrupt dynamic shifts | penalty |
| ruptures BinSeg | Changepoint | Fast approximate | n_bkps |

## Tips

- GAMs with REML estimation are preferred over GCV for smooth parameter estimation; REML is less prone to overfitting
- Always run gam.check() to verify basis dimension is sufficient (k-index should be >= 1.0)
- For comparing conditions, the difference smooth approach (by=is_treated) directly tests trajectory divergence
- Changepoint methods complement GAMs: use GAMs for smooth trends, changepoints for abrupt shifts
- AIC/BIC comparison between linear and GAM fits reveals whether non-linear modeling is warranted

## Related Skills

- temporal-clustering - Group genes after trajectory fitting
- circadian-rhythms - Periodic trajectory models
- differential-expression/timeseries-de - Linear model alternatives for temporal DE
- single-cell/trajectory-inference - Single-cell pseudotime trajectories

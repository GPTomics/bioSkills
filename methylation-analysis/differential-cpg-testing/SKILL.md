---
name: bio-methylation-differential-cpg
description: Per-CpG differential methylation testing from bisulfite sequencing count data or beta-value matrices. Covers beta and M-value computation, coverage filtering, statistical tests (Welch t-test, Mann-Whitney, limma, DSS beta-binomial), multiple testing correction, and effect size calculation. Use when comparing methylation at individual CpG sites between experimental groups from WGBS, RRBS, or targeted bisulfite sequencing.
tool_type: mixed
primary_tool: scipy
---

## Version Compatibility

Reference examples tested with: scipy 1.12+, statsmodels 0.14+, pandas 2.2+, numpy 1.26+, limma 3.58+, DSS 2.50+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Per-CpG Differential Methylation Testing

**"Test individual CpG sites for differential methylation between groups"** -> Compute per-CpG methylation metrics from count data, apply coverage filters, run statistical tests for group differences, correct for multiple testing, and report effect sizes.
- Python: `scipy.stats.ttest_ind()` + `statsmodels.stats.multitest.multipletests()`
- R: `limma::lmFit()` + `eBayes()` on M-values, `DSS::DMLtest()` on counts

## Beta Values and M-Values

**Goal:** Convert raw bisulfite sequencing count data into analyzable methylation metrics.

**Approach:** Compute beta values (methylation proportion) for biological interpretation and M-values (logit-transformed) for statistical testing. Beta values are bounded [0, 1] and heteroscedastic (variance depends on methylation level). M-values are approximately homoscedastic, making them better suited for parametric tests (Du et al. 2010).

```python
import numpy as np
import pandas as pd

# Beta value: proportion of methylated reads
# Range [0, 1], directly interpretable as methylation percentage
beta = meth_counts / total_counts

# M-value: logit transform of beta (log2 scale)
# Range (-inf, +inf), homoscedastic, better for statistical testing
# Add offset to avoid log(0) and division by zero
OFFSET = 1e-6
m_value = np.log2((beta + OFFSET) / (1 - beta + OFFSET))
```

```r
# Beta value from counts
beta <- meth_counts / total_counts

# M-value: logit transform (base 2)
# Add offset to handle beta = 0 or beta = 1
offset <- 1e-6
m_value <- log2((beta + offset) / (1 - beta + offset))

# Convert M-value back to beta for reporting
# beta = 2^M / (2^M + 1)
beta_from_m <- 2^m_value / (2^m_value + 1)
```

When to use each metric:
- **Large n (>10/group), t-test or Mann-Whitney:** Testing on beta values directly is standard and acceptable. Heteroscedasticity has minimal impact with sufficient samples.
- **Small n (3-5/group), limma:** Convert to M-values for statistical testing. Empirical Bayes moderation on M-values produces better-calibrated p-values.
- **Reporting:** Always report delta_beta (difference in mean beta values) for biological interpretation regardless of which metric was used for testing.

## Coverage Filtering

**Goal:** Remove CpGs with unreliable methylation estimates due to insufficient or artificially inflated read coverage.

**Approach:** Apply a minimum coverage threshold to ensure adequate resolution, and an upper percentile filter to remove PCR amplification artifacts.

```python
# Minimum 10x: provides 11 possible methylation levels (0/10 through 10/10)
# and adequate power for statistical testing
MIN_COVERAGE = 10

# Upper 99.9th percentile: removes sites with inflated coverage from
# PCR duplication or mapping artifacts in repetitive regions
upper_threshold = np.percentile(total_counts.values.flatten(), 99.9)

# Require minimum coverage in ALL samples (any low-coverage sample
# makes that CpG unreliable for group comparison)
passes_min = (total_counts >= MIN_COVERAGE).all(axis=1)
passes_max = (total_counts <= upper_threshold).all(axis=1)
filtered = beta[passes_min & passes_max]
```

```r
# Equivalent in R (without methylKit)
min_coverage <- 10
upper_threshold <- quantile(as.matrix(total_counts), 0.999)
keep <- apply(total_counts, 1, function(x) all(x >= min_coverage & x <= upper_threshold))
filtered_beta <- beta[keep, ]
```

### Coverage Thresholds by Assay

| Assay | Minimum | Upper Filter | Rationale |
|-------|---------|--------------|-----------|
| WGBS | 5-10x | 99.9th %ile | Lower thresholds viable when using smoothing methods (BSmooth) |
| RRBS | 10x | 99.9th %ile | Standard minimum; no spatial smoothing for uncovered regions |
| Targeted/amplicon | 30-100x | 99.9th %ile | Deep sequencing expected; higher threshold increases confidence |

## Method Selection

| Scenario | Recommended | Rationale |
|----------|-------------|-----------|
| Large n (>10/group), Python environment | Welch's t-test + BH | Variance estimates reliable at larger sample sizes; fast and simple |
| Large n, non-normal beta distributions | Mann-Whitney U + BH | No distributional assumptions; same power as t-test at large n |
| Small n (3-5/group) | limma on M-values | Empirical Bayes borrows variance across CpGs; adds ~10-20 effective df |
| Count-level modeling needed | DSS beta-binomial | Models both biological variation and sampling noise from read counts |
| Unreplicated (two samples only) | Fisher's exact on counts | Only valid option; cannot estimate biological variance without replicates |
| Mixed tissue samples | limma + cell-type covariates | Cell composition confounds methylation differences (Jaffe & Irizarry 2014) |

Fisher's exact test on counts should be avoided when biological replicates exist. Pooling counts across replicates and testing with Fisher's inflates false positive rates to ~33% because biological variation is ignored.

### Fisher's Exact Test (Python, Unreplicated Only)

```python
from scipy.stats import fisher_exact

# Only for unreplicated designs (1 sample per group)
# table: [[meth_case, unmeth_case], [meth_ctrl, unmeth_ctrl]]
table = [[meth_case, total_case - meth_case],
         [meth_ctrl, total_ctrl - meth_ctrl]]
result = fisher_exact(table, alternative='two-sided')
# result.statistic = odds ratio, result.pvalue = p-value
```

## Welch's t-Test on Beta Values (Python)

**Goal:** Test each CpG for a mean methylation difference between two groups using a parametric test that does not assume equal variances.

**Approach:** Compute beta values per CpG per sample, run a Welch t-test per CpG across groups, collect p-values, and apply Benjamini-Hochberg FDR correction.

```python
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

# Assume: case_betas and ctrl_betas are DataFrames
# Rows = CpGs, columns = samples within each group
pvalues = []
for cpg_idx in range(len(case_betas)):
    case_vals = case_betas.iloc[cpg_idx].values
    ctrl_vals = ctrl_betas.iloc[cpg_idx].values
    # equal_var=False: Welch's t-test (does NOT assume equal variance)
    # scipy defaults to equal_var=True (Student's), which is almost never appropriate
    result = stats.ttest_ind(case_vals, ctrl_vals, equal_var=False, nan_policy='omit')
    pvalues.append(result.pvalue)

pvalues = np.array(pvalues)

# Benjamini-Hochberg FDR correction
# CRITICAL: multipletests defaults to method='hs' (Holm-Sidak), NOT BH
# Must explicitly pass method='fdr_bh'
reject, padj, _, _ = multipletests(pvalues, alpha=0.05, method='fdr_bh')
```

For large datasets, vectorized computation avoids the per-CpG loop:

```python
from scipy.stats import ttest_ind

# Vectorized across all CpGs at once (axis=1 tests across samples)
t_stats, pvalues = ttest_ind(case_betas.values, ctrl_betas.values,
                              axis=1, equal_var=False, nan_policy='omit')

reject, padj, _, _ = multipletests(pvalues, alpha=0.05, method='fdr_bh')
```

## Mann-Whitney U Test (Python)

**Goal:** Non-parametric alternative when beta-value distributions are non-normal or sample sizes are very unequal.

**Approach:** Per-CpG rank-based test with BH correction. No distributional assumptions, but lower power than t-test for small n.

```python
from scipy.stats import mannwhitneyu

pvalues = []
for cpg_idx in range(len(case_betas)):
    case_vals = case_betas.iloc[cpg_idx].dropna().values
    ctrl_vals = ctrl_betas.iloc[cpg_idx].dropna().values
    # Cannot reach p < 0.05 with n < 4 per group (insufficient ranks)
    if len(case_vals) < 4 or len(ctrl_vals) < 4:
        pvalues.append(np.nan)
        continue
    result = mannwhitneyu(case_vals, ctrl_vals, alternative='two-sided')
    pvalues.append(result.pvalue)

pvalues = np.array(pvalues)
valid = ~np.isnan(pvalues)
padj = np.full_like(pvalues, np.nan)
padj[valid] = multipletests(pvalues[valid], method='fdr_bh')[1]
```

## limma on M-Values (R)

**Goal:** Identify differentially methylated CpGs with small sample sizes by borrowing variance information across CpGs via empirical Bayes moderation.

**Approach:** Convert beta values to M-values, fit a linear model per CpG with limma, apply empirical Bayes variance shrinkage, and extract moderated test statistics. This adds ~10-20 effective degrees of freedom even with n=3 per group.

```r
library(limma)

# M-value matrix: rows = CpGs, columns = samples
offset <- 1e-6
m_values <- log2((beta_matrix + offset) / (1 - beta_matrix + offset))

# Design matrix
group <- factor(c(rep('case', n_case), rep('ctrl', n_ctrl)))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# Contrast: case vs control
contrast_matrix <- makeContrasts(case - ctrl, levels = design)

# Fit linear model and apply empirical Bayes
fit <- lmFit(m_values, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
# trend=TRUE: intensity-dependent prior variance (recommended for methylation)
# robust=TRUE: protects against outlier CpGs
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)

# Extract results (all CpGs, BH-adjusted)
results <- topTable(fit2, number = Inf, adjust.method = 'BH', sort.by = 'none')
# Columns: logFC (on M-value scale), AveExpr, t, P.Value, adj.P.Val, B

# Convert logFC (M-value scale) back to approximate delta_beta for reporting
# For CpGs in the middle range (beta 0.2-0.8), delta_beta ~ 0.15 * logFC
# For precise delta_beta, compute from original beta values instead
delta_beta <- rowMeans(beta_matrix[, group == 'case']) -
              rowMeans(beta_matrix[, group == 'ctrl'])
results$delta_beta <- delta_beta
```

## DSS Beta-Binomial Model (R)

**Goal:** Statistically principled count-based differential methylation testing that models both biological variation (Beta distribution across replicates) and sampling noise (Binomial distribution from read counts).

**Approach:** DSS uses a hierarchical beta-binomial model with Bayesian shrinkage for dispersion estimation. It operates directly on methylated/total count data without requiring beta or M-value conversion.

```r
library(DSS)

# Create BSseq object from count matrices
# meth_counts: matrix of methylated read counts (CpGs x samples)
# total_counts: matrix of total read counts (CpGs x samples)
bs_obj <- BSseq(chr = cpg_chr, pos = cpg_pos,
                M = meth_counts, Cov = total_counts,
                sampleNames = sample_names)

# Define groups
group_case <- sample_names[1:n_case]
group_ctrl <- sample_names[(n_case + 1):(n_case + n_ctrl)]

# Test for differential methylation (per CpG)
# smoothing=FALSE: test individual CpGs without spatial smoothing
# smoothing=TRUE: borrows information from nearby CpGs (500bp default window)
dml_test <- DMLtest(bs_obj, group1 = group_case, group2 = group_ctrl, smoothing = FALSE)

# Extract significant DMCs
# p.threshold: FDR-adjusted p-value cutoff
dmcs <- callDML(dml_test, p.threshold = 0.05)
```

## Effect Size Calculation

**Goal:** Quantify the magnitude of methylation difference at each CpG for biological interpretation and downstream filtering.

**Approach:** Compute delta_beta as the difference in mean beta values between groups. Apply effect size thresholds appropriate to the biological context.

```python
delta_beta = case_betas.mean(axis=1).values - ctrl_betas.mean(axis=1).values
```

```r
delta_beta <- rowMeans(beta_matrix[, group == 'case']) -
              rowMeans(beta_matrix[, group == 'ctrl'])
```

### Effect Size Thresholds

| Threshold | delta_beta | Typical Use Case |
|-----------|------------|------------------|
| Lenient | 0.10 | EWAS discovery, environmental exposure studies |
| Standard | 0.20 | General differential methylation analysis |
| methylKit default | 0.25 | methylKit getMethylDiff (reported as 25% difference) |
| Stringent | 0.30 | High-confidence set, cancer vs normal comparisons |

The appropriate threshold depends on biological context. Promoter CpG islands are often bimodal (either <10% or >80% methylated), so even small changes near the transition boundary can be functionally meaningful.

methylKit's `meth.diff` is derived from its logistic regression model and is not identical to `mean(case_beta) - mean(ctrl_beta)`. When comparing results across tools, always recompute delta_beta from the original beta values for consistency.

## Full Python Pipeline

**Goal:** Complete per-CpG differential methylation analysis from a count-data TSV through to a results table with statistics and effect sizes.

**Approach:** Read count data, compute beta values, apply coverage filters, run Welch's t-test per CpG, correct for multiple testing with BH FDR, compute delta_beta, and write results.

```python
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

# Read count data: columns are cpg_id, {sample}_meth, {sample}_total
counts = pd.read_csv('bisulfite_counts.tsv', sep='\t', index_col='cpg_id')

# Identify samples and groups from column naming pattern
case_samples = [c.replace('_meth', '') for c in counts.columns if c.startswith('case') and c.endswith('_meth')]
ctrl_samples = [c.replace('_meth', '') for c in counts.columns if c.startswith('ctrl') and c.endswith('_meth')]

# Compute beta values: meth / total per CpG per sample
case_betas = pd.DataFrame({s: counts[f'{s}_meth'] / counts[f'{s}_total'] for s in case_samples}, index=counts.index)
ctrl_betas = pd.DataFrame({s: counts[f'{s}_meth'] / counts[f'{s}_total'] for s in ctrl_samples}, index=counts.index)

# Coverage filter: require minimum 10x in ALL samples
MIN_COVERAGE = 10  # 11 possible methylation levels, adequate power for testing
total_cols = [c for c in counts.columns if c.endswith('_total')]
passes_coverage = (counts[total_cols] >= MIN_COVERAGE).all(axis=1)
case_betas = case_betas[passes_coverage]
ctrl_betas = ctrl_betas[passes_coverage]

# Per-CpG Welch's t-test (vectorized)
t_stats, pvalues = ttest_ind(case_betas.values, ctrl_betas.values,
                              axis=1, equal_var=False, nan_policy='omit')

# BH FDR correction (default is Holm-Sidak, must specify fdr_bh)
reject, padj, _, _ = multipletests(pvalues, alpha=0.05, method='fdr_bh')

# Effect sizes
mean_case = case_betas.mean(axis=1)
mean_ctrl = ctrl_betas.mean(axis=1)
delta_beta = mean_case - mean_ctrl

# Assemble results
results = pd.DataFrame({
    'cpg_id': case_betas.index,
    'mean_case_beta': mean_case.values,
    'mean_ctrl_beta': mean_ctrl.values,
    'delta_beta': delta_beta.values,
    'pvalue': pvalues,
    'padj': padj,
    'significant': np.where(padj < 0.05, 'TRUE', 'FALSE')
})
results.to_csv('dmc.tsv', sep='\t', index=False)
```

## Common Pitfalls

- **Beta-value heteroscedasticity:** Variance of beta values depends on methylation level -- highest near 0.5, compressed near 0 and 1. This means a 10% difference near beta=0 has different statistical properties than 10% near beta=0.5. For small samples, use limma on M-values instead.

- **scipy.stats.ttest_ind defaults to Student's t-test** (`equal_var=True`). Always pass `equal_var=False` for Welch's t-test. Equal variance between groups is almost never true for methylation data.

- **statsmodels multipletests defaults to Holm-Sidak** (`method='hs'`), not Benjamini-Hochberg. Always pass `method='fdr_bh'` explicitly.

- **Fisher's exact test with replicates inflates false positives.** Pooling methylated/unmethylated counts across biological replicates and running Fisher's test ignores biological variation, producing false positive rates of ~33%. Fisher's exact is only appropriate for unreplicated designs (single sample per group).

- **methylKit meth.diff is not mean(case_beta) - mean(ctrl_beta).** methylKit uses logistic regression with overdispersion correction. Its `meth.diff` column reflects the model-estimated difference, which can diverge from the simple mean difference, particularly at extreme methylation levels.

- **Bonferroni correction is too conservative for genome-wide methylation.** With 450K-850K CpGs (arrays) or 28M+ CpGs (WGBS), Bonferroni thresholds leave almost no significant CpGs. BH FDR is the standard approach.

- **Cell-type composition confounds mixed-tissue samples.** In blood or bulk tissue, observed methylation differences may reflect cell-type proportion changes rather than within-cell methylation changes (Jaffe & Irizarry 2014). Include cell-type estimates as covariates or use deconvolution methods.

## Related Skills

- methylkit-analysis - Chi-squared testing with methylKit in R
- methylation-calling - Generate per-CpG count data from Bismark BAM files
- dmr-detection - Region-level testing after per-CpG analysis
- differential-expression/deseq2-basics - Analogous empirical Bayes concepts for count data
- proteomics/differential-abundance - Similar Welch t-test + BH workflow for proteomics

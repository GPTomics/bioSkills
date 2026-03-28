---
name: bio-proteomics-differential-abundance
description: Statistical testing for differentially abundant proteins between conditions. Covers preprocessing (log2 transformation, normalization), limma and DEqMS workflows with empirical Bayes moderation, fold change shrinkage for accurate effect size estimation, and Python alternatives. Use when identifying proteins with significant abundance changes between experimental groups.
tool_type: mixed
primary_tool: limma
---

## Version Compatibility

Reference examples tested with: limma 3.58+, DEqMS 1.20+, ashr 2.2+, proDA 1.20+, numpy 1.26+, pandas 2.2+, scipy 1.12+, statsmodels 0.14+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Differential Protein Abundance

**"Find differentially abundant proteins between my conditions"** -> Perform statistical testing on protein intensities to identify significant abundance changes.
- R: `limma::eBayes()` for empirical Bayes moderated t-tests (preferred for small n)
- R: `DEqMS::spectraCounteBayes()` when PSM/peptide count metadata is available
- R: `proDA::test_diff()` when missing values are extensive (label-free)
- Python: `scipy.stats.ttest_ind(equal_var=False)` with `statsmodels` BH correction

## Preprocessing Pipeline

Raw mass spectrometry intensities require log2 transformation and normalization before statistical testing.

### Log2 Transformation

**Goal:** Convert right-skewed raw intensities to approximately normal distributions with stabilized variance.

**Approach:** Apply log2 to all intensity values. Replace zeros (undetected values) with NaN before transformation to avoid -inf.

```python
log2_data = np.log2(intensities.replace(0, np.nan))
```

```r
log2_matrix <- log2(intensity_matrix)
log2_matrix[!is.finite(log2_matrix)] <- NA
```

### Normalization

**Goal:** Remove systematic technical biases (sample loading, instrument drift) so that observed differences reflect biology.

**Approach:** Choose a normalization method based on data characteristics. All methods assume the majority of proteins are not differentially abundant.

**Median normalization** -- subtract per-sample median so all samples share a common center:

```python
sample_medians = log2_data.median(axis=0)
global_median = sample_medians.median()
normalized = log2_data - sample_medians + global_median
```

```r
normalized <- normalizeBetweenArrays(log2_matrix, method = 'scale')
```

| Method | When to use | R function |
|--------|-------------|------------|
| Median centering | Default for most analyses; robust to missing values | Manual or `normalizeBetweenArrays(method='scale')` |
| Cyclic loess | Unbalanced DE (more up- than down-regulated) | `normalizeBetweenArrays(method='cyclicloess')` |
| VSN | Heteroscedastic data; operates on raw intensities (skip log2) | `vsn::justvsn(raw_matrix)` |
| Quantile | TMT with complete data; avoid with many missing values | `normalizeBetweenArrays(method='quantile')` |

## Method Selection

| Scenario | Recommended | Rationale |
|----------|-------------|-----------|
| Small n (3-5 per group), protein-level data | limma | Borrows variance across proteins via empirical Bayes; adds ~10-20 effective df |
| PSM/peptide count metadata available | DEqMS | Weights variance by quantification depth per protein |
| Label-free with many missing values (>20%) | proDA | Models abundance-dependent dropout; no imputation needed |
| Large n (>10 per group), Python-only environment | Welch's t-test + BH | Variance estimates reliable at larger sample sizes |
| Complex designs (nested, multiple comparisons) | MSstats | Feature-level mixed models; handles technical replicates |

With small sample sizes (n=3-5), simple t-tests have only 4-6 degrees of freedom for variance estimation, making per-protein variances extremely noisy. Some proteins get artificially low variance (false positives), others artificially high (false negatives). limma's empirical Bayes shrinks each variance toward the global trend, dramatically improving calibration.

## limma Workflow (R)

**Goal:** Identify differentially abundant proteins using moderated statistics that borrow information across all proteins.

**Approach:** Build a linear model from the design matrix, fit contrasts, apply empirical Bayes moderation with intensity-dependent variance trend and robust fitting, then extract BH-corrected results.

```r
library(limma)

design <- model.matrix(~0 + condition, data = sample_info)
colnames(design) <- levels(factor(sample_info$condition))

fit <- lmFit(protein_matrix, design)
contrast_matrix <- makeContrasts(Treatment - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)

results <- topTable(fit2, coef = 1, number = Inf, adjust.method = 'BH')
```

`trend=TRUE` allows the prior variance to depend on mean intensity. `robust=TRUE` protects against hyper-variable outlier proteins. Results columns: `logFC`, `AveExpr`, `t`, `P.Value`, `adj.P.Val`, `B`.

For batch effects, include batch in the design matrix (do NOT use `removeBatchEffect` before testing -- that function is for visualization only):

```r
design <- model.matrix(~0 + condition + batch, data = sample_info)
```

## DEqMS Workflow (R)

**Goal:** Improve upon limma by accounting for the relationship between quantification depth and variance -- proteins quantified by more PSMs/peptides have more precise abundance estimates.

**Approach:** Run the standard limma pipeline, attach PSM/peptide counts, then apply DEqMS's count-aware empirical Bayes that fits a variance-vs-count regression.

```r
library(DEqMS)

# Standard limma pipeline through eBayes (see above), then:
fit2$count <- psm_count_per_protein[rownames(fit2$coefficients)]
fit3 <- spectraCounteBayes(fit2)

results <- outputResult(fit3, coef_col = 1)
```

Results include limma columns plus DEqMS-specific: `sca.t`, `sca.P.Value`, `sca.adj.pval`.

## proDA Workflow (R)

```r
library(proDA)

fit <- proDA(protein_matrix, design = ~condition, col_data = sample_info,
             reference_level = 'Control')
results <- test_diff(fit, conditionTreatment - conditionControl)
```

Results columns: `name`, `pval`, `adj_pval`, `diff` (log2FC), `t_statistic`, `se`.

## Python Workflow

**Goal:** Perform the full differential abundance pipeline in Python: preprocessing, statistical testing, and multiple testing correction.

**Approach:** Log2-transform and median-normalize raw intensities, run per-protein Welch's t-tests, and apply Benjamini-Hochberg correction.

```python
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

def preprocess(intensities):
    log2_data = np.log2(intensities.replace(0, np.nan))
    sample_medians = log2_data.median(axis=0)
    global_median = sample_medians.median()
    return log2_data - sample_medians + global_median

def differential_abundance(normalized, case_cols, ctrl_cols):
    results = []
    for protein in normalized.index:
        case = normalized.loc[protein, case_cols].dropna()
        ctrl = normalized.loc[protein, ctrl_cols].dropna()
        if len(case) >= 2 and len(ctrl) >= 2:
            log2fc = case.mean() - ctrl.mean()
            _, pval = stats.ttest_ind(case, ctrl, equal_var=False)
            results.append({'protein': protein, 'log2fc': log2fc, 'pvalue': pval})

    df = pd.DataFrame(results)
    df['padj'] = multipletests(df['pvalue'], method='fdr_bh')[1]
    return df
```

**Key details:**
- `equal_var=False` selects Welch's t-test (`scipy` defaults to Student's with `equal_var=True`)
- `multipletests` defaults to Holm-Sidak -- always pass `method='fdr_bh'` explicitly

## Fold Change Reporting

Raw fold changes from the linear model or t-test are the best unbiased point estimates of the true effect, but they are noisy -- proteins with no real abundance change still show small nonzero estimates from measurement noise. How to handle this depends on the downstream use case.

### When to report raw fold changes

Report the unmodified log2 fold change from the statistical test when:
- Running **GSEA or pathway analysis** that ranks all proteins by effect size (these methods rely on the full continuous distribution, including small non-significant effects)
- Performing **meta-analysis** across studies (raw FCs with standard errors are the correct input)
- The downstream consumer needs an **unbiased estimate** with associated uncertainty (report FC + SE or confidence interval)

### When to apply fold change shrinkage

Apply shrinkage when the goal is to **recover which proteins truly changed and by how much** -- i.e., effect size accuracy matters more than preserving the full distribution. ashr (R) fits a mixture prior with a point mass at zero and estimates posterior means, smoothly shrinking uncertain effects toward zero while preserving well-supported ones:

```r
library(ashr)

se <- sqrt(fit2$s2.post) * fit2$stdev.unscaled[, 1]
shrunk <- ash(fit2$coefficients[, 1], se, mixcompdist = 'normal')

shrunken_fc <- shrunk$result$PosteriorMean
lfsr <- shrunk$result$lfsr
```

ashr is preferred over hard-thresholding (zeroing FCs at a p-value cutoff) because it shrinks smoothly based on per-protein uncertainty rather than applying an arbitrary step function at padj = 0.05. Hard zeroing discards information and creates artificial discontinuities -- a protein at padj 0.049 keeps its full FC while one at 0.051 is set to zero.

In Python without ashr, there is no mature equivalent. If effect size accuracy is critical, use R with ashr. For Python-only environments, report raw fold changes with adjusted p-values and let the downstream analysis handle thresholding.

### Minimum fold change testing

To test whether fold changes exceed a biologically meaningful threshold (rather than just differ from zero), use `treat()` + `topTreat()` instead of post-hoc FC filtering with `topTable(lfc=...)`, which can inflate FDR:

```r
fit2 <- treat(fit2, lfc = log2(1.2))
results <- topTreat(fit2, coef = 1, number = Inf)
```

## Visualization

```r
library(ggplot2)

ggplot(results, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(color = significant), alpha = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
    geom_vline(xintercept = c(-1, 1), linetype = 'dashed') +
    scale_color_manual(values = c('grey60', 'firebrick')) +
    theme_minimal() + labs(x = 'Log2 Fold Change', y = '-Log10 Adjusted P-value')
```

## Common Pitfalls

- **Not log-transforming raw intensities** -- parametric tests assume approximately normal distributions; raw intensities are right-skewed with mean-dependent variance
- **Using Student's t-test** -- `scipy.stats.ttest_ind` defaults to `equal_var=True`; always set `equal_var=False` (Welch's) since treatment can affect both mean and variance
- **Quantile normalization with missing values** -- introduces artifacts in label-free data; use median centering or cyclic loess instead
- **`removeBatchEffect()` before testing** -- this function is for visualization only; include batch as a covariate in the design matrix for statistical testing
- **Post-hoc FC filtering via `topTable(lfc=...)`** -- can inflate FDR; use `treat()` + `topTreat()` for minimum-effect-size testing
- **Ignoring fold change uncertainty** -- raw FCs are noisy point estimates; consider ashr shrinkage when effect size accuracy matters, and always report adjusted p-values or confidence intervals alongside fold changes so downstream analyses can weight accordingly

## Related Skills

- quantification - Protein-level abundance estimation and normalization before testing
- proteomics-qc - Quality control and batch effect assessment
- differential-expression/deseq2-basics - Analogous empirical Bayes concepts for RNA-seq
- data-visualization/specialized-omics-plots - Volcano plots, MA plots, heatmaps

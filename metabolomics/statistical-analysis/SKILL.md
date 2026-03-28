---
name: bio-metabolomics-statistical-analysis
description: Statistical analysis for metabolomics data. Covers preprocessing (log2 transformation, normalization), limma moderated testing with empirical Bayes, Welch's t-tests with BH correction, fold change estimation, and multivariate methods (PCA, PLS-DA, OPLS-DA). Use when identifying differentially abundant metabolites or building classification models.
tool_type: mixed
primary_tool: limma
---

## Version Compatibility

Reference examples tested with: limma 3.58+, ashr 2.2+, scipy 1.12+, statsmodels 0.14+, numpy 1.26+, pandas 2.1+, mixOmics 6.28+, ggplot2 3.5+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Metabolomics Statistical Analysis

## Processing Pipeline

Standard pipeline: zero/missing handling -> log2 transformation -> normalization -> statistical testing. Fold change and statistical tests operate on log2-transformed, normalized but unscaled data. Pareto or auto-scaling distorts fold change magnitudes -- apply only for multivariate methods (PCA, PLS-DA).

### Zero and Missing Value Handling

Metabolomics zeros fall into two categories requiring different treatment:

| Type | Meaning | Detection | Action |
|------|---------|-----------|--------|
| TPMV (technical) | Below detection limit | Random missingness within detected features | Impute: half-minimum per feature or KNN |
| BPMV (biological) | Metabolite truly absent | Structured: all zeros in one group | Leave as-is or use two-part test |

Distinguishing them matters: imputing BPMVs introduces false signal, while ignoring TPMVs loses power. If zeros are concentrated in one experimental group, they are likely biological. If scattered randomly, they are likely technical.

```python
log2_data = np.log2(intensities.replace(0, np.nan))
```

```r
log2_matrix <- log2(intensity_matrix)
log2_matrix[!is.finite(log2_matrix)] <- NA
```

For pseudocount approach: `log2(x + 1)` avoids NaN but compresses low-intensity fold changes. Prefer half-minimum imputation when the zero rate is low (<10%).

### Normalization

**Goal:** Remove systematic technical biases (injection volume, instrument drift, ionization efficiency) so observed differences reflect biology.

**Approach:** Choose based on data characteristics. All methods assume the majority of features are not differentially abundant. Unlike proteomics where median centering is the default, metabolomics uses PQN because it is more robust to dominant high-abundance features common in MS-based metabolomics.

| Method | When to use | Metabolomics note |
|--------|-------------|-------------------|
| PQN | **Default for untargeted LC-MS metabolomics** | More robust than TIC/median to dominant features |
| QC-RSC (LOESS) | Multi-batch LC-MS with pooled QC samples | Gold standard for batch correction; metabolomics-specific |
| VSN | High zero rate or heteroscedastic data | Handles zeros via arcsinh; replaces separate log2 step |
| TIC | Quick exploration; NMR data | Distorted by dominant features; avoid for LC-MS |
| Cyclic loess | Asymmetric DE (more up- than down-regulated) | Robust to assumption violations |
| None | IS-corrected data; single-batch balanced design | Check PCA on raw log2 data first |

PQN normalization: compute a reference spectrum (median across samples), divide each feature by the reference, take the median of quotients as the per-sample normalization factor.

```python
reference = log2_data.median(axis=1)
quotients = log2_data.div(reference, axis=0)
norm_factors = quotients.median(axis=0)
normalized = log2_data.div(norm_factors, axis=1)
```

```r
reference <- apply(log2_matrix, 1, median, na.rm = TRUE)
quotients <- sweep(log2_matrix, 1, reference, '/')
norm_factors <- apply(quotients, 2, median, na.rm = TRUE)
normalized <- sweep(log2_matrix, 2, norm_factors, '/')
```

**When normalization hurts:** Drug studies and extreme phenotypes can cause global metabolic shifts affecting >50% of features, violating the "majority stable" assumption. This is more common in metabolomics than proteomics because metabolic networks are tightly coupled. Check PCA on raw log2 data first -- if dominant variation corresponds to the biological factor AND global shifts are expected, skip or attenuate normalization.

## Method Selection

| Scenario | Recommended | Rationale |
|----------|-------------|-----------|
| Small n (3-10/group), default | limma `eBayes(trend=TRUE, robust=TRUE)` | Borrows variance across features; adds ~10-20 effective df |
| Large n (>10/group), Python-only | Welch's t-test + BH | Per-feature variance reliable; limma converges to ordinary t-test |
| Non-normal after log transform | Wilcoxon rank-sum | Cannot reach p < 0.05 with n < 4/group |
| Zero-inflated (many BPMVs) | Two-part test | Separately tests presence/absence and abundance |
| Paired/repeated measures | Paired t-test or limma with blocking | `duplicateCorrelation()` in limma for repeated measures |
| 3+ groups | Welch's ANOVA or limma F-test | Post-hoc: Games-Howell or Tukey HSD |

Unlike proteomics, metabolomics has no analog to PSM/peptide counts, so DEqMS is not applicable. limma is the sole empirical Bayes option for metabolite-level testing.

`trend=TRUE` is particularly important for metabolomics: features span 3-4 orders of magnitude in intensity and include chemically diverse compounds (polar metabolites, lipids, nucleotides), creating a stronger mean-variance relationship than typical proteomics data. `robust=TRUE` protects against hyper-variable outlier features common in metabolomics panels that mix chemical classes.

## limma Workflow (R)

**Goal:** Identify differentially abundant metabolites using moderated statistics that borrow information across all features.

**Approach:** Build a linear model, apply empirical Bayes moderation with intensity-dependent variance trend, extract BH-corrected results.

```r
library(limma)

design <- model.matrix(~0 + condition, data = sample_info)
colnames(design) <- levels(factor(sample_info$condition))

fit <- lmFit(normalized_matrix, design)
contrast_matrix <- makeContrasts(Treatment - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)

results <- topTable(fit2, coef = 1, number = Inf, adjust.method = 'BH')
```

Results columns: `logFC`, `AveExpr`, `t`, `P.Value`, `adj.P.Val`, `B`. Note: adjusted p-values are in `adj.P.Val` (not `FDR` or `padj`). For batch effects, include batch in the design matrix (do NOT use `removeBatchEffect()` before testing -- that function is for visualization only):

```r
design <- model.matrix(~0 + condition + batch, data = sample_info)
```

## Welch's t-test Workflow (Python)

**Goal:** Perform the full differential abundance pipeline in Python when R/limma is unavailable.

**Approach:** Log2-transform raw intensities, run per-feature Welch's t-tests, apply BH correction. Add PQN normalization when systematic biases are suspected (see Processing Pipeline).

```python
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

intensities = pd.read_csv('feature_table.tsv', sep='\t', index_col=0)
metadata = pd.read_csv('sample_info.tsv', sep='\t')
case_samples = metadata.loc[metadata['group'] == 'case', 'sample_id'].tolist()
ctrl_samples = metadata.loc[metadata['group'] == 'control', 'sample_id'].tolist()

log2_data = np.log2(intensities.replace(0, np.nan))

pvalues, log2fcs = [], []
for feature in log2_data.index:
    case_vals = log2_data.loc[feature, case_samples].dropna().values.astype(float)
    ctrl_vals = log2_data.loc[feature, ctrl_samples].dropna().values.astype(float)
    if len(case_vals) >= 2 and len(ctrl_vals) >= 2:
        _, pval = ttest_ind(case_vals, ctrl_vals, equal_var=False)  # Welch's -- scipy defaults to Student's
        pvalues.append(pval)
        log2fcs.append(case_vals.mean() - ctrl_vals.mean())
    else:
        pvalues.append(np.nan)
        log2fcs.append(np.nan)

results = pd.DataFrame({'feature_id': log2_data.index, 'log2fc': log2fcs, 'pvalue': pvalues})
results = results.dropna(subset=['pvalue'])
_, results['padj'], _, _ = multipletests(results['pvalue'], method='fdr_bh')  # default is Holm-Sidak, not BH
results['significant'] = results['padj'] < 0.05
```

**Key details:**
- `equal_var=False` selects Welch's t-test (scipy defaults to Student's with `equal_var=True`)
- `multipletests` defaults to Holm-Sidak (FWER) -- always pass `method='fdr_bh'` explicitly for BH FDR
- Fold change is computed as difference of group means in log2-space (geometric mean ratio), not `log2(mean_ratio)` which uses arithmetic means

## Fold Change Calculation

| Scenario | Test | Notes |
|----------|------|-------|
| 2 groups, n >= 5/group | Welch's t-test | Always prefer over Student's; unequal variance is the norm |
| 2 groups, non-normal after log | Mann-Whitney U | Cannot reach p < 0.05 with n < 4/group |
| 2 groups, n < 5/group | limma moderated t | `eBayes(trend=TRUE)` borrows variance across features |
| Paired samples | Paired t-test | Pre/post, matched case-control |
| 3+ groups | Welch's ANOVA | Post-hoc: Games-Howell or Dunn's test |

**Approach:** Compute the difference of group means on log2-transformed data, which equals log2(geometric_mean_case / geometric_mean_control).

```python
log2fc = log2_data.loc[:, case_samples].mean(axis=1) - log2_data.loc[:, ctrl_samples].mean(axis=1)
```

```r
log2fc <- rowMeans(normalized[, case_samples]) - rowMeans(normalized[, ctrl_samples])
```

Always compute fold change on log2-transformed, unscaled data. The naive alternative -- `log2(mean(case) / mean(control))` -- uses arithmetic means and can reverse fold change direction when group variances differ. The difference-of-log-means approach uses geometric means, consistent with limma and DESeq2.

**MetaboAnalyst note:** MetaboAnalyst internally uses `log2(arithmetic_mean_ratio)` for fold change calculation. This diverges from the geometric mean approach and can produce different results, particularly for features with heterogeneous variance. The geometric mean method (difference of log means) is statistically more appropriate for log-normally distributed metabolomics data.

### Fold Change Reporting

Raw fold changes are noisy but unbiased. How to handle them depends on the downstream use:

- **Pathway analysis (QEA/GSEA)**: Use raw fold changes for all features. MetaboAnalyst QEA ranks by effect size and relies on the full continuous distribution. Do not zero or threshold FCs before enrichment analysis.
- **Effect size recovery** (which metabolites truly changed and by how much): Apply ashr shrinkage (R) for posterior mean estimates. Relevant for biomarker panels or clinical reporting.
- **Meta-analysis across studies**: Use raw fold changes with standard errors as input. Shrinkage is study-specific and should not be applied before pooling.

### Fold Change Shrinkage (ashr)

ashr fits a mixture prior with a point mass at zero and estimates posterior means, smoothly shrinking uncertain effects toward zero while preserving well-supported ones:

```r
library(ashr)

se <- sqrt(fit2$s2.post) * fit2$stdev.unscaled[, 1]
shrunk <- ash(fit2$coefficients[, 1], se, mixcompdist = 'normal')

shrunken_fc <- shrunk$result$PosteriorMean
lfsr <- shrunk$result$lfsr
```

ashr is preferred over hard-thresholding (zeroing FCs at a p-value cutoff) because it shrinks smoothly based on per-feature uncertainty rather than applying an arbitrary step function at padj = 0.05. In Python, there is no mature equivalent -- use R with ashr when effect size accuracy is critical.

### Minimum Fold Change Testing

To test whether fold changes exceed a biologically meaningful threshold rather than just differ from zero, use `treat()` + `topTreat()` instead of post-hoc FC filtering with `topTable(lfc=...)`, which can inflate FDR:

```r
fit2 <- treat(fit2, lfc = log2(1.5))
results <- topTreat(fit2, coef = 1, number = Inf)
```

## Common Pitfalls

- **Not log-transforming raw intensities** -- parametric tests assume approximately normal distributions; raw MS intensities are right-skewed with mean-dependent variance
- **Using Student's t-test** -- `scipy.stats.ttest_ind` defaults to `equal_var=True`; always set `equal_var=False` (Welch's) since treatment can affect both mean and variance
- **statsmodels default is not BH** -- `multipletests` defaults to Holm-Sidak (FWER); always pass `method='fdr_bh'` for metabolomics FDR correction
- **`removeBatchEffect()` before testing** -- this function is for visualization only; include batch as a covariate in the design matrix for statistical testing
- **Post-hoc FC filtering** -- `topTable(lfc=...)` inflates FDR; use `treat()` + `topTreat()` for minimum-effect-size testing
- **Normalizing when majority-stable fails** -- drug studies causing global metabolic shifts violate PQN/TIC/median assumptions; check PCA on raw log2 data first
- **Imputing biological zeros** -- replacing BPMVs (metabolite truly absent in one group) with half-minimum creates false signal; use two-part tests instead
- **MetaboAnalyst defaults** -- uses Student's t-test (not Welch's) and arithmetic mean FC (not geometric mean) by default; override or be aware when comparing results

## Volcano Plot

**Goal:** Visualize both statistical significance and effect size for all features in a single plot.

**Approach:** Plot log2 fold change (x-axis) vs -log10 p-value (y-axis), highlighting features passing both thresholds.

```python
import matplotlib.pyplot as plt

results['sig_and_fc'] = (results['padj'] < 0.05) & (results['log2fc'].abs() > 1)
colors = ['red' if s else 'gray' for s in results['sig_and_fc']]

plt.figure(figsize=(8, 6))
plt.scatter(results['log2fc'], -np.log10(results['pvalue']), c=colors, alpha=0.6, s=20)
plt.axhline(-np.log10(0.05), linestyle='--', color='gray')
plt.axvline(-1, linestyle='--', color='gray')
plt.axvline(1, linestyle='--', color='gray')
plt.xlabel('Log2 Fold Change')
plt.ylabel('-log10(p-value)')
plt.savefig('volcano_plot.png', dpi=150, bbox_inches='tight')
```

```r
library(ggplot2)

results$sig_and_fc <- results$adj.P.Val < 0.05 & abs(results$logFC) > 1

ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), color = sig_and_fc)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c('gray60', 'firebrick')) +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
    geom_vline(xintercept = c(-1, 1), linetype = 'dashed') +
    labs(x = 'Log2 Fold Change', y = '-Log10 Adjusted P-value') +
    theme_bw()
```

## PCA

**Goal:** Explore overall sample variation and detect outliers or batch effects in an unsupervised manner.

**Approach:** Perform PCA on the feature matrix and plot the first two principal components colored by experimental group.

```r
library(pcaMethods)

pca_result <- pca(data, nPcs = 5, method = 'ppca')

scores <- as.data.frame(scores(pca_result))
scores$group <- groups

ggplot(scores, aes(x = PC1, y = PC2, color = group)) +
    geom_point(size = 3) +
    stat_ellipse(level = 0.95) +
    labs(x = paste0('PC1 (', round(pca_result@R2[1] * 100, 1), '%)'),
         y = paste0('PC2 (', round(pca_result@R2[2] * 100, 1), '%)')) +
    theme_bw()

loadings <- as.data.frame(loadings(pca_result))
top_pc1 <- loadings[order(abs(loadings$PC1), decreasing = TRUE)[1:20], ]
```

## PLS-DA

**Goal:** Build a supervised classification model that maximizes separation between experimental groups.

**Approach:** Fit a PLS-DA model with cross-validation to determine optimal components, then extract VIP scores to rank discriminatory features.

```r
library(mixOmics)

plsda_result <- plsda(as.matrix(data), groups, ncomp = 3)

perf_plsda <- perf(plsda_result, validation = 'Mfold', folds = 5, nrepeat = 50)
plot(perf_plsda, col = color.mixo(5:7))

ncomp_opt <- perf_plsda$choice.ncomp['BER', 'centroids.dist']

final_plsda <- plsda(as.matrix(data), groups, ncomp = ncomp_opt)
plotIndiv(final_plsda, group = groups, ellipse = TRUE, legend = TRUE)

vip <- vip(final_plsda)
top_vip <- sort(vip[, ncomp_opt], decreasing = TRUE)[1:20]
```

## sPLS-DA (Sparse)

**Goal:** Perform feature selection simultaneously with classification to identify a minimal discriminatory feature set.

**Approach:** Tune the number of features to retain per component via cross-validation, then fit a sparse PLS-DA model.

```r
tune_splsda <- tune.splsda(as.matrix(data), groups, ncomp = 3,
                            validation = 'Mfold', folds = 5, nrepeat = 50,
                            test.keepX = c(5, 10, 20, 50, 100))

optimal_keepX <- tune_splsda$choice.keepX

splsda_result <- splsda(as.matrix(data), groups, ncomp = ncomp_opt, keepX = optimal_keepX)

selected_features <- selectVar(splsda_result, comp = 1)$name
```

## OPLS-DA

**Goal:** Separate group-predictive variation from orthogonal (within-group) variation for cleaner class separation.

**Approach:** Fit an OPLS-DA model using ropls, then use the S-plot to identify features with high predictive weight and correlation.

```r
library(ropls)

oplsda <- opls(data, groups, predI = 1, orthoI = NA)
plot(oplsda, typeVc = 'x-score')
plot(oplsda, typeVc = 'x-loading')

vip_scores <- getVipVn(oplsda)
top_vip <- sort(vip_scores, decreasing = TRUE)[1:20]
```

## Random Forest

**Goal:** Rank features by importance using a non-linear ensemble classifier.

**Approach:** Train a Random Forest on the feature matrix, extract MeanDecreaseAccuracy to identify discriminatory features.

```r
library(randomForest)

rf_model <- randomForest(x = data, y = groups, importance = TRUE, ntree = 500)
importance <- importance(rf_model)
top_features <- rownames(importance)[order(importance[, 'MeanDecreaseAccuracy'], decreasing = TRUE)[1:20]]
varImpPlot(rf_model, n.var = 20)
```

## ROC Analysis

**Goal:** Evaluate the diagnostic performance of candidate biomarker metabolites.

**Approach:** Generate ROC curves and compute AUC for individual features using pROC.

```r
library(pROC)

top_feature <- 'feature_123'
roc_result <- roc(groups, data[, top_feature])
plot(roc_result, main = paste('AUC =', round(auc(roc_result), 3)))

biomarkers <- c('feature_1', 'feature_2', 'feature_3')
for (feat in biomarkers) {
    roc_i <- roc(groups, data[, feat])
    cat(feat, ': AUC =', round(auc(roc_i), 3), '\n')
}
```

## Heatmap

**Goal:** Visualize abundance patterns of top differential features across all samples.

**Approach:** Select top significant features, create an annotated heatmap with hierarchical clustering.

```r
library(pheatmap)

top_features <- rownames(sig_features)[1:50]
data_top <- data[, top_features]
annotation_row <- data.frame(Group = groups)
rownames(annotation_row) <- rownames(data)

pheatmap(t(data_top), annotation_col = annotation_row,
         scale = 'row', clustering_method = 'ward.D2',
         filename = 'heatmap.png', width = 10, height = 12)
```

## Related Skills

- normalization-qc - Data preparation and batch correction
- pathway-mapping - Functional interpretation of differential metabolites
- multi-omics-integration/mixomics-analysis - Advanced multivariate methods
- proteomics/differential-abundance - Analogous empirical Bayes concepts for proteomics
- data-visualization/specialized-omics-plots - Volcano plots, MA plots, heatmaps

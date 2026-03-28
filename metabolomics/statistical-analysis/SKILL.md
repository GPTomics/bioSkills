---
name: bio-metabolomics-statistical-analysis
description: Statistical analysis for metabolomics data. Covers univariate testing, multivariate methods (PCA, PLS-DA), and biomarker discovery. Use when identifying differentially abundant metabolites or building classification models.
tool_type: mixed
primary_tool: mixOmics
---

## Version Compatibility

Reference examples tested with: scipy 1.12+, statsmodels 0.14+, numpy 1.26+, pandas 2.1+, R stats (base), mixOmics 6.28+, ggplot2 3.5+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Metabolomics Statistical Analysis

## Processing Order

Standard pipeline: missing value handling → sample normalization (TIC/PQN) → log2 transformation → statistical testing. Fold change and t-tests operate on log2-transformed, unscaled data. Pareto or auto-scaling distorts fold change magnitudes — apply only for multivariate methods (PCA, PLS-DA).

If raw data contains zeros, impute before log transformation: half-minimum per feature (`min(nonzero) / 2`) for below-LOD values, or use a pseudocount (`log2(x + 1)`).

## Univariate Analysis

**Goal:** Identify differentially abundant metabolites between experimental groups using per-feature statistical tests on log2-transformed data.

**Approach:** Log2-transform raw intensities, apply Welch's t-tests per feature, then correct for multiple testing with Benjamini-Hochberg FDR.

**"Find differentially abundant metabolites between my groups"** → Log2-transform intensities, apply per-feature Welch's t-tests, correct with BH FDR.
- Python: `scipy.stats.ttest_ind(equal_var=False)` + `statsmodels.stats.multitest.multipletests()`
- R: `t.test()` (Welch's by default) + `p.adjust(method='BH')`

```python
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

intensities = pd.read_csv('feature_table.csv', index_col=0)
metadata = pd.read_csv('sample_info.csv')
treat_samples = metadata.loc[metadata['group'] == 'treatment', 'sample_id'].tolist()
ctrl_samples = metadata.loc[metadata['group'] == 'control', 'sample_id'].tolist()

log2_data = np.log2(intensities)

pvalues, log2fcs = [], []
for feature in log2_data.index:
    treat_vals = log2_data.loc[feature, treat_samples].values.astype(float)
    ctrl_vals = log2_data.loc[feature, ctrl_samples].values.astype(float)
    _, pval = ttest_ind(treat_vals, ctrl_vals, equal_var=False)  # Welch's -- scipy defaults to Student's
    pvalues.append(pval)
    log2fcs.append(treat_vals.mean() - ctrl_vals.mean())

results = pd.DataFrame({'feature_id': log2_data.index, 'log2fc': log2fcs, 'pvalue': pvalues})
_, results['padj'], _, _ = multipletests(results['pvalue'], method='fdr_bh')
results['significant'] = results['padj'] < 0.05
```

```r
data <- read.csv('feature_table.csv', row.names = 1)
groups <- factor(read.csv('sample_info.csv')$group)

data_log2 <- log2(data)

ttest_results <- apply(data_log2, 2, function(x) {
    treat_vals <- x[groups == 'treatment']
    ctrl_vals <- x[groups == 'control']
    test <- t.test(treat_vals, ctrl_vals)  # Welch's by default (var.equal=FALSE)
    c(pvalue = test$p.value, log2fc = mean(treat_vals) - mean(ctrl_vals))
})
ttest_results <- as.data.frame(t(ttest_results))
ttest_results$padj <- p.adjust(ttest_results$pvalue, method = 'BH')
sig_features <- ttest_results[ttest_results$padj < 0.05, ]
```

### Test Selection

| Scenario | Test | Notes |
|----------|------|-------|
| 2 groups, n >= 5/group | Welch's t-test | Always prefer over Student's; unequal variance is the norm |
| 2 groups, non-normal after log | Mann-Whitney U | Cannot reach p < 0.05 with n < 4/group |
| 2 groups, n < 5/group | limma moderated t | `eBayes(trend=TRUE)` borrows variance across features |
| Paired samples | Paired t-test | Pre/post, matched case-control |
| 3+ groups | Welch's ANOVA | Post-hoc: Games-Howell or Dunn's test |

## Fold Change Calculation

**Goal:** Quantify the magnitude and direction of abundance changes between groups.

**Approach:** Compute the difference of group means on log2-transformed data, which equals log2(geometric_mean_treatment / geometric_mean_control).

```python
log2_data = np.log2(intensities)
log2fc = log2_data.loc[:, treat_samples].mean(axis=1) - log2_data.loc[:, ctrl_samples].mean(axis=1)
```

```r
data_log2 <- log2(data)
log2fc <- colMeans(data_log2[groups == 'treatment', ]) - colMeans(data_log2[groups == 'control', ])
```

Always compute fold change on log2-transformed, unscaled data. The naive alternative — `log2(mean(case) / mean(control))` — uses arithmetic means and can reverse fold change direction when group variances differ. The difference-of-log-means approach uses geometric means, consistent with limma and DESeq2.

## Volcano Plot

**Goal:** Visualize both statistical significance and effect size for all features in a single plot.

**Approach:** Plot log2 fold change (x-axis) vs -log10 p-value (y-axis), highlighting features passing both thresholds.

```python
import matplotlib.pyplot as plt

results['significant'] = (results['padj'] < 0.05) & (results['log2fc'].abs() > 1)
colors = ['red' if s else 'gray' for s in results['significant']]

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

ttest_results$significant <- ttest_results$padj < 0.05 & abs(ttest_results$log2fc) > 1

ggplot(ttest_results, aes(x = log2fc, y = -log10(pvalue), color = significant)) +
    geom_point(alpha = 0.6) +
    scale_color_manual(values = c('gray', 'red')) +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
    geom_vline(xintercept = c(-1, 1), linetype = 'dashed') +
    labs(x = 'Log2 Fold Change', y = '-log10(p-value)') +
    theme_bw()

ggsave('volcano_plot.png', width = 8, height = 6)
```

## PCA

**Goal:** Explore overall sample variation and detect outliers or batch effects in an unsupervised manner.

**Approach:** Perform PCA on the feature matrix and plot the first two principal components colored by experimental group.

```r
library(pcaMethods)

# PCA
pca_result <- pca(data, nPcs = 5, method = 'ppca')

# Scores plot
scores <- as.data.frame(scores(pca_result))
scores$group <- groups

ggplot(scores, aes(x = PC1, y = PC2, color = group)) +
    geom_point(size = 3) +
    stat_ellipse(level = 0.95) +
    labs(x = paste0('PC1 (', round(pca_result@R2[1] * 100, 1), '%)'),
         y = paste0('PC2 (', round(pca_result@R2[2] * 100, 1), '%)')) +
    theme_bw()

ggsave('pca_plot.png', width = 8, height = 6)

# Loadings
loadings <- as.data.frame(loadings(pca_result))
top_pc1 <- loadings[order(abs(loadings$PC1), decreasing = TRUE)[1:20], ]
```

## PLS-DA

**Goal:** Build a supervised classification model that maximizes separation between experimental groups.

**Approach:** Fit a PLS-DA model with cross-validation to determine optimal components, then extract VIP scores to rank discriminatory features.

```r
library(mixOmics)

# PLS-DA
plsda_result <- plsda(as.matrix(data), groups, ncomp = 3)

# Cross-validation
perf_plsda <- perf(plsda_result, validation = 'Mfold', folds = 5, nrepeat = 50)
plot(perf_plsda, col = color.mixo(5:7))

# Optimal components
ncomp_opt <- perf_plsda$choice.ncomp['BER', 'centroids.dist']
cat('Optimal components:', ncomp_opt, '\n')

# Final model
final_plsda <- plsda(as.matrix(data), groups, ncomp = ncomp_opt)

# Plot
plotIndiv(final_plsda, group = groups, ellipse = TRUE, legend = TRUE)

# VIP scores
vip <- vip(final_plsda)
top_vip <- sort(vip[, ncomp_opt], decreasing = TRUE)[1:20]
print(top_vip)
```

## sPLS-DA (Sparse)

**Goal:** Perform feature selection simultaneously with classification to identify a minimal discriminatory feature set.

**Approach:** Tune the number of features to retain per component via cross-validation, then fit a sparse PLS-DA model.

```r
# Tune number of features to select
tune_splsda <- tune.splsda(as.matrix(data), groups, ncomp = 3,
                            validation = 'Mfold', folds = 5, nrepeat = 50,
                            test.keepX = c(5, 10, 20, 50, 100))

optimal_keepX <- tune_splsda$choice.keepX

# Final sparse model
splsda_result <- splsda(as.matrix(data), groups, ncomp = ncomp_opt, keepX = optimal_keepX)

# Selected features
selected_features <- selectVar(splsda_result, comp = 1)$name
cat('Selected features (comp 1):', length(selected_features), '\n')
```

## OPLS-DA (Orthogonal PLS-DA)

**Goal:** Separate group-predictive variation from orthogonal (within-group) variation for cleaner class separation.

**Approach:** Fit an OPLS-DA model using ropls, then use the S-plot to identify features with high predictive weight and correlation.

```r
library(ropls)

# OPLS-DA
oplsda <- opls(data, groups, predI = 1, orthoI = NA)

# Summary
print(oplsda)

# Scores plot
plot(oplsda, typeVc = 'x-score')

# S-plot (loadings vs correlation)
plot(oplsda, typeVc = 'x-loading')

# VIP
vip_scores <- getVipVn(oplsda)
top_vip <- sort(vip_scores, decreasing = TRUE)[1:20]
```

## Random Forest

**Goal:** Rank features by importance using a non-linear ensemble classifier.

**Approach:** Train a Random Forest on the feature matrix, then extract MeanDecreaseAccuracy to identify the most discriminatory features.

```r
library(randomForest)

# Random Forest classification
rf_model <- randomForest(x = data, y = groups, importance = TRUE, ntree = 500)

# Importance
importance <- importance(rf_model)
top_features <- rownames(importance)[order(importance[, 'MeanDecreaseAccuracy'], decreasing = TRUE)[1:20]]

# Plot importance
varImpPlot(rf_model, n.var = 20)
```

## ROC Analysis

**Goal:** Evaluate the diagnostic performance of candidate biomarker metabolites.

**Approach:** Generate ROC curves and compute AUC for individual features using pROC.

```r
library(pROC)

# ROC for top biomarker
top_feature <- 'feature_123'  # Replace with actual feature name
roc_result <- roc(groups, data[, top_feature])

# Plot
plot(roc_result, main = paste('AUC =', round(auc(roc_result), 3)))

# Multiple features
biomarkers <- c('feature_1', 'feature_2', 'feature_3')
for (feat in biomarkers) {
    roc_i <- roc(groups, data[, feat])
    cat(feat, ': AUC =', round(auc(roc_i), 3), '\n')
}
```

## Heatmap

**Goal:** Visualize abundance patterns of top differential features across all samples.

**Approach:** Select top significant features, create an annotated heatmap with hierarchical clustering using pheatmap.

```r
library(pheatmap)

# Top differential features
top_features <- rownames(sig_features)[1:50]
data_top <- data[, top_features]

# Annotation
annotation_row <- data.frame(Group = groups)
rownames(annotation_row) <- rownames(data)

pheatmap(t(data_top), annotation_col = annotation_row,
         scale = 'row', clustering_method = 'ward.D2',
         filename = 'heatmap.png', width = 10, height = 12)
```

## Related Skills

- normalization-qc - Data preparation
- pathway-mapping - Functional interpretation
- multi-omics-integration/mixomics-analysis - Advanced multivariate methods

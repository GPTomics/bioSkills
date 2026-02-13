---
name: bio-workflows-timecourse-pipeline
description: End-to-end time-course analysis from expression matrix to temporal patterns and enrichment. Covers temporal DE, Mfuzz soft clustering, optional rhythm detection, GAM trajectory fitting, and per-cluster pathway enrichment. Use when analyzing bulk time-series expression experiments from any omics platform.
tool_type: mixed
primary_tool: Mfuzz
workflow: true
depends_on:
  - differential-expression/timeseries-de
  - temporal-genomics/temporal-clustering
  - temporal-genomics/circadian-rhythms
  - temporal-genomics/trajectory-modeling
  - pathway-analysis/go-enrichment
qc_checkpoints:
  - after_de: "Significant temporal genes >100 at FDR <0.05; model fit residuals reasonable"
  - after_clustering: "Membership >0.5 for soft clustering; no empty clusters; validated by gap statistic or silhouette (typical range 4-20 clusters depending on complexity)"
  - after_enrichment: "At least 3 clusters with significant GO terms at FDR <0.05"
---

## Version Compatibility

Reference examples tested with: DESeq2 1.42+, clusterProfiler 4.10+, limma 3.58+, numpy 1.26+, pandas 2.2+, scanpy 1.10+, scipy 1.12+, statsmodels 0.14+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Time-Course Analysis Pipeline

**"Analyze my time-course expression data end-to-end"** → Orchestrate temporal differential expression, Mfuzz soft clustering, optional circadian rhythm detection, GAM trajectory fitting, changepoint detection, and per-cluster pathway enrichment.

Complete workflow from expression matrix through temporal differential expression, soft clustering,
optional rhythm detection, trajectory fitting, and per-cluster pathway enrichment.

## Pipeline Overview

```
Expression matrix + time metadata
    |
    v
[1. Temporal DE] ---------> limma splines / DESeq2 LRT
    |
    v
[2. Filter] --------------> Significant temporal genes (FDR <0.05)
    |
    v
[3. Mfuzz Clustering] ----> Soft clustering of expression profiles
    |                            |
    |                            +---> QC: membership >0.5, no empty clusters
    |
    +--- Circadian design? ---> [4a. Rhythm Detection] (MetaCycle / CosinorPy)
    |                               |
    |                               v
    |                           Rhythmic genes + period/phase estimates
    |
    v
[4b. GAM Trajectory] -----> mgcv GAM fitting for top clusters
    |
    v
[5. Pathway Enrichment] --> clusterProfiler per-cluster GO/KEGG
    |
    v
Temporal gene modules + enriched pathways + trajectory plots
```

## Step 1: Temporal Differential Expression

### R (limma splines)

```r
library(limma)
library(splines)

expr <- as.matrix(read.csv('counts_normalized.csv', row.names = 1))
meta <- read.csv('metadata.csv')

time_points <- meta$time
design <- model.matrix(~ ns(time_points, df = 3))

fit <- lmFit(expr, design)
fit <- eBayes(fit)

# Test all spline coefficients jointly for temporal significance
temporal_results <- topTable(fit, coef = 2:ncol(design), number = Inf, sort.by = 'F')
# topTable already returns adj.P.Val (BH-corrected); use it directly
```

### R (DESeq2 LRT)

```r
library(DESeq2)

counts <- as.matrix(read.csv('raw_counts.csv', row.names = 1))
meta <- read.csv('metadata.csv')
meta$time <- factor(meta$time)

dds <- DESeqDataSetFromMatrix(counts, colData = meta, design = ~ time)
dds <- DESeq(dds, test = 'LRT', reduced = ~ 1)
res <- results(dds)
```

### Python (statsmodels)

```python
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests
from patsy import dmatrix
from scipy import stats

expr = pd.read_csv('counts_normalized.csv', index_col=0)
meta = pd.read_csv('metadata.csv')

spline_basis = dmatrix('bs(time, df=3)', data=meta, return_type='dataframe')
design_full = np.column_stack([np.ones(len(meta)), spline_basis.values])
design_reduced = np.ones((len(meta), 1))

pvals = []
for gene in expr.index:
    y = expr.loc[gene].values
    ss_full = np.sum((y - design_full @ np.linalg.lstsq(design_full, y, rcond=None)[0]) ** 2)
    ss_red = np.sum((y - design_reduced @ np.linalg.lstsq(design_reduced, y, rcond=None)[0]) ** 2)
    df_diff = design_full.shape[1] - design_reduced.shape[1]
    df_resid = len(y) - design_full.shape[1]
    f_stat = ((ss_red - ss_full) / df_diff) / (ss_full / df_resid)
    pvals.append(1 - stats.f.cdf(f_stat, df_diff, df_resid))

_, fdr, _, _ = multipletests(pvals, method='fdr_bh')
temporal_genes = expr.index[fdr < 0.05].tolist()
```

### QC Checkpoint: Temporal DE

```r
# Gate 1: Sufficient temporal genes detected
sig_genes <- temporal_results[temporal_results$adj.P.Val < 0.05, ]
n_sig <- nrow(sig_genes)
message(sprintf('Significant temporal genes: %d', n_sig))
if (n_sig < 100) message('WARNING: Few temporal genes. Check time point spacing or consider relaxing FDR.')
if (n_sig > 10000) message('WARNING: Many temporal genes. Consider stricter FDR or inspect batch effects.')

# Gate 2: Residual distribution check
residuals <- residuals(fit, expr)
message(sprintf('Residual mean: %.4f, SD: %.4f', mean(residuals), sd(residuals)))
```

## Step 2: Filter Significant Genes

```r
# FDR <0.05: Standard threshold for temporal DE
# More permissive (0.1) acceptable for exploratory clustering
sig_genes <- rownames(temporal_results[temporal_results$adj.P.Val < 0.05, ])
expr_sig <- expr[sig_genes, ]
message(sprintf('Genes passing FDR <0.05: %d', length(sig_genes)))
```

## Step 3: Mfuzz Soft Clustering

```r
library(Mfuzz)

eset <- ExpressionSet(assayData = as.matrix(expr_sig))

# Standardize expression profiles (mean=0, sd=1 per gene)
eset <- standardise(eset)

# Estimate fuzzifier m
# mestimate() calculates optimal m from data geometry; typical range 1.5-2.5
m <- mestimate(eset)
message(sprintf('Estimated fuzzifier m = %.2f', m))

# Cluster count: start with sqrt(n_genes/2), refine with gap statistic
# Typical range 4-20 depending on temporal complexity
n_clusters <- 8
cl <- mfuzz(eset, c = n_clusters, m = m)

# Filter low-membership genes
# Membership >0.5: gene clearly belongs to one cluster
# Lower to 0.3 for exploratory analysis with more overlap
core_genes <- acore(eset, cl, min.acore = 0.5)
```

### Python Alternative (tslearn)

```python
from tslearn.clustering import TimeSeriesKMeans

# Row-wise z-scoring: normalize each gene across its timepoints (not per-timepoint)
expr_scaled = (expr_sig.values - expr_sig.values.mean(axis=1, keepdims=True)) / expr_sig.values.std(axis=1, keepdims=True)

# n_clusters: 4-20 depending on complexity; evaluate with silhouette score
model = TimeSeriesKMeans(n_clusters=8, metric='softdtw', metric_params={'gamma': 0.01},
                         max_iter=50, random_state=42)
labels = model.fit_predict(expr_scaled.reshape(expr_scaled.shape[0], expr_scaled.shape[1], 1))
```

### QC Checkpoint: Clustering

```r
# Gate 1: No empty clusters
cluster_sizes <- table(cl$cluster)
message('Cluster sizes:')
print(cluster_sizes)
if (any(cluster_sizes == 0)) message('WARNING: Empty clusters found. Reduce n_clusters.')

# Gate 2: Membership quality
for (i in seq_along(core_genes)) {
    n_core <- nrow(core_genes[[i]])
    message(sprintf('Cluster %d: %d core genes (membership >0.5)', i, n_core))
}

# Gate 3: Silhouette score (optional validation)
library(cluster)
hard_labels <- cl$cluster
sil <- silhouette(hard_labels, dist(exprs(eset)))
message(sprintf('Mean silhouette: %.3f', mean(sil[, 3])))
```

## Step 4a: Rhythm Detection (Optional - Circadian Designs)

Only applicable when sampling covers 24h+ cycles with sufficient resolution (every 2-4h).

### R (MetaCycle)

```r
library(MetaCycle)

# Expects genes as rows, time points as columns
# Column names must be numeric time values (hours)
expr_for_meta <- expr_sig
colnames(expr_for_meta) <- meta$time_hours

write.csv(expr_for_meta, 'expr_for_metacycle.csv')

# Period range 20-28h: standard circadian search window
# Adjust for ultradian (4-12h) or infradian (>28h) rhythms
meta2d('expr_for_metacycle.csv', filestyle = 'csv',
       minper = 20, maxper = 28,
       timepoints = sort(unique(meta$time_hours)),
       outdir = 'metacycle_results')
```

### Python (CosinorPy)

```python
from cosinorpy import file_parser, cosinor

# fit_group expects long-format DataFrame with columns 'x' (time), 'y' (expression), 'test' (gene name)
# Reshape expression matrix to long format before passing
# period=24: standard circadian; adjust for other periodicities
results = cosinor.fit_group(expr_long, period=24, n_components=1)
rhythmic = results[results['p'] < 0.05]
```

## Step 4b: GAM Trajectory Fitting

### R (mgcv)

```r
library(mgcv)

cluster_trajectories <- list()
for (cl_id in 1:n_clusters) {
    cl_genes <- names(cl$cluster[cl$cluster == cl_id])
    mean_profile <- colMeans(expr_sig[cl_genes, ])

    df_gam <- data.frame(time = meta$time, expr = mean_profile)

    # k: basis dimension; k=5 sufficient for most time courses
    # Increase to k=10 for >20 time points; decrease to k=3 for <6 time points
    gam_fit <- gam(expr ~ s(time, k = 5), data = df_gam)

    cluster_trajectories[[cl_id]] <- list(
        fit = gam_fit,
        r_squared = summary(gam_fit)$r.sq,
        edf = summary(gam_fit)$edf
    )
    message(sprintf('Cluster %d: R^2 = %.3f, EDF = %.2f', cl_id,
                    summary(gam_fit)$r.sq, summary(gam_fit)$edf))
}
```

### Python (pygam)

```python
from pygam import LinearGAM, s
import numpy as np

for cl_id in range(n_clusters):
    cl_mask = labels == cl_id
    mean_profile = expr_scaled[cl_mask].mean(axis=0)

    # n_splines=5: sufficient for most time courses
    gam = LinearGAM(s(0, n_splines=5)).fit(meta['time'].values.reshape(-1, 1), mean_profile)
    print(f'Cluster {cl_id}: GCV = {gam.statistics_["GCV"]:.4f}')
```

## Step 5: Per-Cluster Pathway Enrichment

### R (clusterProfiler)

```r
library(clusterProfiler)
library(org.Hs.eg.db)

enrichment_results <- list()
for (i in seq_along(core_genes)) {
    genes <- core_genes[[i]]$NAME

    # Convert symbols to Entrez IDs
    entrez <- bitr(genes, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)

    # GO Biological Process enrichment
    # pvalueCutoff 0.05: standard; qvalueCutoff 0.05 for FDR control
    ego <- enrichGO(gene = entrez$ENTREZID, OrgDb = org.Hs.eg.db,
                    ont = 'BP', pAdjustMethod = 'BH',
                    pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                    readable = TRUE)

    enrichment_results[[i]] <- ego
    message(sprintf('Cluster %d: %d significant GO terms', i, nrow(as.data.frame(ego))))
}
```

### Python (gseapy)

```python
import gseapy as gp

for cl_id in range(n_clusters):
    cl_genes = [g for g, l in zip(expr_sig.index, labels) if l == cl_id]

    enr = gp.enrichr(gene_list=cl_genes, gene_sets='GO_Biological_Process_2023',
                     organism='human', outdir=f'enrichr_cluster_{cl_id}')
    sig_terms = enr.results[enr.results['Adjusted P-value'] < 0.05]
    print(f'Cluster {cl_id}: {len(sig_terms)} significant GO terms')
```

### QC Checkpoint: Enrichment

```r
# Gate: At least 3 clusters should have significant GO terms
clusters_with_terms <- sum(sapply(enrichment_results, function(x) nrow(as.data.frame(x)) > 0))
message(sprintf('Clusters with significant GO terms: %d / %d', clusters_with_terms, length(enrichment_results)))
if (clusters_with_terms < 3) {
    message('WARNING: Few clusters enriched. Check gene ID conversion or relax thresholds.')
}
```

## Parameter Recommendations

| Step | Parameter | Recommendation |
|------|-----------|----------------|
| Temporal DE | Spline df | 3 (default); increase to 4-5 for >10 time points |
| Temporal DE | FDR | 0.05 (standard); 0.1 for exploratory clustering |
| Mfuzz | fuzzifier m | Use mestimate(); typical range 1.5-2.5 |
| Mfuzz | n_clusters | 4-20; start with sqrt(n_genes/2), refine with gap statistic |
| Mfuzz | min membership | 0.5 (core genes); 0.3 (exploratory) |
| MetaCycle | period range | 20-28h (circadian); adjust for other periodicities |
| GAM | k (basis dim) | 5 (default); 3 for <6 time points; 10 for >20 |
| clusterProfiler | pvalueCutoff | 0.05 (standard); 0.1 (permissive) |

## Troubleshooting

| Issue | Likely Cause | Solution |
|-------|--------------|----------|
| < 100 temporal genes | Insufficient replicates or noisy data | Add replicates; use DESeq2 LRT instead of limma |
| Empty Mfuzz clusters | Too many clusters | Reduce n_clusters; check gap statistic |
| All genes in one cluster | Fuzzifier too low or too few clusters | Increase m or n_clusters |
| No rhythmic genes | Non-circadian design or low power | Verify 24h+ sampling; increase resolution |
| GAM overfitting | k too high for time points | Set k = min(n_timepoints - 1, 5) |
| Few enriched clusters | Gene ID conversion failure | Check species; verify Entrez ID mapping |
| Low membership scores | High expression noise | Increase fuzzifier m; apply stricter gene filtering |

## Related Skills

- differential-expression/timeseries-de - Temporal DE methods
- temporal-genomics/temporal-clustering - Cluster analysis details
- temporal-genomics/circadian-rhythms - Rhythm detection details
- temporal-genomics/trajectory-modeling - GAM fitting details
- pathway-analysis/go-enrichment - Enrichment analysis details
- temporal-genomics/periodicity-detection - Unknown-period discovery

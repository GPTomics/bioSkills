---
name: bio-temporal-genomics-temporal-clustering
description: Clusters genes by temporal expression profile shape using Mfuzz soft clustering, TCseq, and DEGreport degPatterns. Groups co-regulated genes into shared trajectory patterns via fuzzy c-means or hierarchical approaches. Use when categorizing temporally dynamic genes into response groups or identifying co-expression modules across time points. Requires temporally variable genes identified first (see differential-expression/timeseries-de).
tool_type: mixed
primary_tool: Mfuzz
---

## Version Compatibility

Reference examples tested with: numpy 1.26+, scanpy 1.10+, scikit-learn 1.4+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Temporal Gene Clustering

**"Group my time-course genes by expression pattern shape"** → Cluster temporally variable genes into co-expression modules by trajectory shape using fuzzy c-means (Mfuzz), hierarchical methods, or DTW-based approaches, revealing coordinated response patterns.
- R: `Mfuzz::mfuzz()` for soft (fuzzy) temporal clustering
- Python: `sklearn.cluster.KMeans` on z-scored time profiles for hard clustering

Groups genes with similar temporal expression dynamics into clusters, revealing shared regulatory programs and coordinated response patterns across time-course experiments.

## Core Workflow

1. Select temporally variable genes (pre-filtered by DE or variance)
2. Standardize expression profiles (z-score across timepoints)
3. Choose clustering method and number of clusters
4. Assign genes to clusters (hard or soft membership)
5. Validate clusters and run functional enrichment per cluster

## Mfuzz (R/Bioconductor)

**Goal:** Group temporally variable genes into co-expression clusters by trajectory shape using fuzzy c-means, revealing shared regulatory programs.

**Approach:** Create an ExpressionSet from the time-series matrix, filter low-variance genes, standardize profiles, estimate the fuzzifier parameter, then run fuzzy c-means to assign soft cluster memberships.

Soft (fuzzy) c-means clustering assigns genes membership scores across all clusters, capturing genes with ambiguous temporal behavior.

### Setup and Preprocessing

```r
library(Mfuzz)
library(Biobase)

# Rows = genes, columns = timepoints (mean across replicates)
expr_mat <- as.matrix(read.csv('temporal_expression.csv', row.names = 1))

# Create ExpressionSet
eset <- ExpressionSet(assayData = expr_mat)

# filter.std removes genes with near-zero variance across timepoints
# min.std=0.5: removes flat genes; adjust based on data spread
eset <- filter.std(eset, min.std = 0.5)

# Standardize each gene to mean=0, sd=1 across timepoints
eset <- standardise(eset)
```

### Fuzzifier Estimation and Clustering

```r
# mestimate(): data-driven fuzzifier estimate based on gene count and dimensions
# Prevents clusters from being too crisp (m close to 1) or too fuzzy (m >> 2)
m <- mestimate(eset)
cat(sprintf('Estimated fuzzifier: %.2f\n', m))

# c=8: typical starting point for 6-12 timepoints; refine with cluster validity indices
cl <- mfuzz(eset, c = 8, m = m)

# Membership filtering: genes with membership < 0.5 in all clusters are ambiguous
# 0.5 threshold: standard cutoff; genes below this are equidistant from multiple centroids
core_genes <- acore(eset, cl, min.acore = 0.5)
```

### Visualization

```r
# Temporal profile plot with membership-based color intensity
mfuzz.plot2(eset, cl, mfrow = c(2, 4), time.labels = colnames(expr_mat),
            centre = TRUE, x11 = FALSE)

# Cluster overlap plot shows similarity between cluster centroids
overlap.plot(cl, over = overlap(cl), thres = 0.05)
```

### Cluster Number Selection

```r
# Evaluate multiple k values; pick where cluster validity stabilizes
# Range 4-20: typical for temporal data; fewer for simple designs, more for dense sampling
validity_scores <- numeric()
for (k in 4:20) {
    cl_k <- mfuzz(eset, c = k, m = m)
    # Minimum centroid distance: should not collapse below threshold
    centroids <- cl_k$centers
    dists <- as.matrix(dist(centroids))
    diag(dists) <- Inf
    validity_scores <- c(validity_scores, min(dists))
}
plot(4:20, validity_scores, type = 'b', xlab = 'Number of clusters', ylab = 'Min centroid distance')
```

## TCseq (R/Bioconductor)

Temporal clustering with fuzzy c-means and k-means on time-course sequencing data.

```r
library(TCseq)

# timeclust with fuzzy c-means
# algo='cm': fuzzy c-means; captures soft membership like Mfuzz
# k=6: number of clusters; test range and evaluate with silhouette
tc <- timeclust(expr_mat, algo = 'cm', k = 6, standardize = TRUE)

# Cluster assignment plot
timeclustplot(tc, value = 'z-score', cols = 3)

# k-means alternative for hard clustering
tc_km <- timeclust(expr_mat, algo = 'km', k = 6, standardize = TRUE)
```

## DEGreport degPatterns (R)

Automatic cluster number selection and publication-ready plots.

```r
library(DEGreport)

# degPatterns automatically selects optimal cluster count via hierarchical clustering
# time: factor defining timepoint order
# col: column in metadata for coloring (e.g., condition)
# minc=15: minimum genes per cluster to retain; prevents singleton clusters
patterns <- degPatterns(expr_mat, metadata = sample_info,
                        time = 'timepoint', col = 'condition', minc = 15)

# Access cluster assignments
cluster_df <- patterns$df

# Plot individual clusters
degPlotCluster(patterns$normalized, time = 'timepoint', color = 'condition')
```

## tslearn (Python)

Time-series clustering with Dynamic Time Warping (DTW) distance.

```python
import numpy as np
from tslearn.clustering import TimeSeriesKMeans
from tslearn.preprocessing import TimeSeriesScalerMeanVariance
from tslearn.utils import to_time_series_dataset
from sklearn.metrics import silhouette_score

# expr_mat: numpy array of shape (n_genes, n_timepoints)
expr_scaled = TimeSeriesScalerMeanVariance().fit_transform(expr_mat[:, :, np.newaxis])

# DTW metric: handles phase-shifted profiles better than Euclidean
# Soft-DTW (metric='softdtw') is differentiable and faster for large datasets
# n_clusters=8: starting point; evaluate with silhouette
model = TimeSeriesKMeans(n_clusters=8, metric='dtw', max_iter=50, random_state=42)
labels = model.fit_predict(expr_scaled)

# sklearn silhouette_score does not support DTW; precompute distance matrix
from tslearn.metrics import cdist_dtw
dist_matrix = cdist_dtw(expr_scaled)
sil = silhouette_score(dist_matrix, labels, metric='precomputed')
```

### Cluster Number Selection with Silhouette

```python
# Test k from 3-15; pick k with highest silhouette score
# 3-15 range: fewer than 3 is too coarse; more than 15 rarely adds biological meaning
sil_scores = []
for k in range(3, 16):
    model = TimeSeriesKMeans(n_clusters=k, metric='softdtw', max_iter=30, random_state=42)
    labels = model.fit_predict(expr_scaled)
    # Euclidean silhouette as computational shortcut; DTW silhouette is O(n^2 * T^2)
    sil_scores.append(silhouette_score(expr_scaled.squeeze(), labels, metric='euclidean'))
```

## Method Comparison

| Method | Clustering Type | Distance | Best For |
|--------|----------------|----------|----------|
| Mfuzz | Soft (fuzzy c-means) | Euclidean | Standard temporal profiling |
| TCseq | Soft or hard | Euclidean | RNA-seq time courses |
| DEGreport | Hierarchical | Correlation | Automatic k selection |
| tslearn | Hard (k-means) | DTW/soft-DTW | Phase-shifted profiles |

## Tips

- Always standardize (z-score) before clustering; otherwise, highly expressed genes dominate
- Soft clustering (Mfuzz) is preferred when genes may participate in multiple temporal programs
- DTW-based clustering captures time-shifted patterns but is computationally expensive for >5000 genes
- Run functional enrichment (GO/GSEA) per cluster to interpret biological meaning
- Membership threshold of 0.5 for Mfuzz filters ~30-50% of genes as ambiguous; adjust if too stringent

## Related Skills

- circadian-rhythms - Rhythm-specific clustering by phase
- trajectory-modeling - Continuous trajectory fitting before clustering
- differential-expression/timeseries-de - Upstream temporal DE for gene selection
- pathway-analysis/go-enrichment - Per-cluster functional enrichment

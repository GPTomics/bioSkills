---
name: bio-temporal-genomics-temporal-grn
description: Infers dynamic gene regulatory networks from bulk time-series expression data using Granger causality (statsmodels), dynGENIE3 (Extra-Trees on ODE-derived expression derivatives), and dynamic Bayesian networks (bnlearn). Identifies time-delayed regulatory relationships and tracks network rewiring across conditions. Use when inferring causal regulatory relationships from bulk temporal expression data or detecting TF influence propagation over time. Not for static co-expression networks (see gene-regulatory-networks/coexpression-networks).
tool_type: mixed
primary_tool: statsmodels
---

# Temporal Gene Regulatory Network Inference

Infers directed, time-delayed regulatory relationships from bulk time-series expression data. Captures how transcription factor activity propagates through gene regulatory networks over time.

## Core Workflow

1. Select candidate regulators (transcription factors) and target genes
2. Prepare lagged expression matrices from time-series data
3. Apply temporal inference method (Granger, dynGENIE3, or DBN)
4. Filter significant regulatory edges by statistical threshold
5. Compare networks across conditions or time windows

## Granger Causality (Python/statsmodels)

Tests whether past values of gene X improve prediction of gene Y beyond past values of Y alone. Significant Granger causality suggests X temporally influences Y.

### Pairwise Granger Test

```python
import numpy as np
import pandas as pd
from statsmodels.tsa.stattools import grangercausalitytests
from statsmodels.stats.multitest import multipletests

# expr_df: genes x timepoints DataFrame (must be temporally ordered)
tf_genes = ['TF1', 'TF2', 'TF3']
target_genes = ['geneA', 'geneB', 'geneC']

# maxlag=2: test lags 1 and 2; for 4h sampling, lag 1 = 4h delay, lag 2 = 8h
# Higher maxlag captures longer delays but requires more timepoints (need n > 3*maxlag)
maxlag = 2

results = []
for tf in tf_genes:
    for target in target_genes:
        if tf == target:
            continue
        pair_data = pd.DataFrame({'target': expr_df.loc[target], 'tf': expr_df.loc[tf]}).values
        gc_result = grangercausalitytests(pair_data, maxlag=maxlag, verbose=False)
        # Use minimum p-value across tested lags
        min_p = min(gc_result[lag][0]['ssr_ftest'][1] for lag in range(1, maxlag + 1))
        best_lag = min(range(1, maxlag + 1), key=lambda l: gc_result[l][0]['ssr_ftest'][1])
        results.append({'tf': tf, 'target': target, 'p_value': min_p, 'best_lag': best_lag})

results_df = pd.DataFrame(results)
```

### Multiple Testing Correction

```python
# BH FDR correction: standard for genome-wide regulatory edge testing
reject, qvals, _, _ = multipletests(results_df['p_value'], method='fdr_bh')
results_df['q_value'] = qvals

# q < 0.05: standard FDR threshold for regulatory edge significance
significant_edges = results_df[results_df['q_value'] < 0.05].copy()
significant_edges = significant_edges.sort_values('q_value')
```

### Build Adjacency Matrix

```python
all_genes = list(set(tf_genes + target_genes))
adj_matrix = pd.DataFrame(0.0, index=all_genes, columns=all_genes)

for _, row in significant_edges.iterrows():
    # -log10(q): edge weight; higher = more significant regulatory relationship
    adj_matrix.loc[row['tf'], row['target']] = -np.log10(row['q_value'])
```

### Stationarity Check

```python
from statsmodels.tsa.stattools import adfuller

# Granger causality assumes stationarity; apply uniform first differencing
# Uniform differencing avoids mixing differenced and non-differenced genes
expr_df = expr_df.diff(axis=1).iloc[:, 1:]
```

## dynGENIE3 (R)

Extension of GENIE3 for time-series data. Estimates expression derivatives from consecutive timepoints, then uses Extra-Trees regression to predict derivatives from current expression of candidate regulators.

### Basic dynGENIE3

```r
library(dynGENIE3)

# TS.data: list of expression matrices (one per time series/replicate)
# Each matrix: genes x timepoints
# time.points: corresponding list of time vectors
expr_list <- list(as.matrix(expr_series1), as.matrix(expr_series2))
time_list <- list(c(0, 4, 8, 12, 24, 48), c(0, 4, 8, 12, 24, 48))

# regulators: indices of known TF genes (restricts source nodes)
# Using TF indices improves both accuracy and speed
tf_indices <- which(rownames(expr_list[[1]]) %in% tf_names)

res <- dynGENIE3(
    TS.data = expr_list,
    time.points = time_list,
    regulators = tf_indices
)

# getLinkList: extracts ranked TF-target pairs with importance scores
# threshold=1000: top 1000 edges; adjust based on number of TFs
# Typical: 10-50 edges per TF for focused networks
link_list <- get.link.list(res$weight.matrix, report.max = 1000)
head(link_list)
```

### With Known Regulators

```r
# When TF list is available, restrict regulators for cleaner results
# AnimalTFDB or PlantTFDB provide species-specific TF lists
tf_names <- readLines('tf_list.txt')
tf_idx <- which(rownames(expr_list[[1]]) %in% tf_names)

res_tf <- dynGENIE3(
    TS.data = expr_list,
    time.points = time_list,
    regulators = tf_idx
)
```

### Multiple Time Series

```r
# Multiple biological replicates improve ODE derivative estimation
# Each replicate provides independent derivative samples
expr_list <- list(replicate1_mat, replicate2_mat, replicate3_mat)
time_list <- list(timepoints, timepoints, timepoints)

res_multi <- dynGENIE3(
    TS.data = expr_list,
    time.points = time_list,
    regulators = tf_idx
)
```

## Dynamic Bayesian Networks (R/bnlearn)

Models probabilistic dependencies between genes across time slices.

### Structure Learning

```r
library(bnlearn)

# Encode temporal lags: create columns for t and t-1
# lag=1: one-step dependency; each gene at time t depends on genes at t-1
lagged_df <- data.frame()
for (t in 2:ncol(expr_mat)) {
    row <- c(expr_mat[, t], expr_mat[, t - 1])
    lagged_df <- rbind(lagged_df, row)
}
colnames(lagged_df) <- c(
    paste0(rownames(expr_mat), '_t'),
    paste0(rownames(expr_mat), '_t1')
)

# Hill-climbing structure learning with BIC score
# score='bic-g': BIC for Gaussian networks; penalizes complexity
dag <- hc(lagged_df, score = 'bic-g')

# Extract temporal edges (t-1 -> t only; discard within-slice edges)
edges <- arcs(dag)
temporal_edges <- edges[grepl('_t1$', edges[, 'from']) & grepl('_t$', edges[, 'to']), ]
```

### Bootstrap for Edge Confidence

```r
# boot.strength: resamples data and learns structure repeatedly
# R=200: number of bootstrap replicates; 200 is standard for moderate datasets
boot_res <- boot.strength(lagged_df, algorithm = 'hc',
                          algorithm.args = list(score = 'bic-g'), R = 200)

# strength >= 0.7: edge appears in >70% of bootstrap networks; reasonably confident
# direction >= 0.5: edge direction is consistent in >50% of appearances
confident_edges <- boot_res[boot_res$strength >= 0.7 & boot_res$direction >= 0.5, ]
```

## Network Comparison Across Conditions

### Jaccard Similarity of Edge Sets

```python
def edge_jaccard(edges_a, edges_b):
    set_a = set(map(tuple, edges_a))
    set_b = set(map(tuple, edges_b))
    intersection = len(set_a & set_b)
    union = len(set_a | set_b)
    if union == 0:
        return 0.0
    return intersection / union

# Compare networks from two conditions
jaccard = edge_jaccard(condition1_edges, condition2_edges)
# Jaccard < 0.3: substantial network rewiring between conditions
# Jaccard > 0.7: networks are largely conserved
```

### Differential Edge Detection

```python
# Identify edges unique to each condition
edges_1 = set(map(tuple, condition1_edges))
edges_2 = set(map(tuple, condition2_edges))

gained_edges = edges_2 - edges_1
lost_edges = edges_1 - edges_2
conserved_edges = edges_1 & edges_2
```

## Method Comparison

| Method | Approach | Strengths | Limitations |
|--------|----------|-----------|-------------|
| Granger causality | VAR model comparison | Statistically principled, pairwise | Assumes stationarity, linear |
| dynGENIE3 | ODE + tree regression | Non-linear, multi-regulator | No p-values, ranking-based |
| DBN (bnlearn) | Probabilistic graphical model | Multivariate, captures conditional dependencies | Computationally expensive |

## Parameter Guide

| Parameter | Typical Value | Rationale |
|-----------|---------------|-----------|
| Granger maxlag | 1-3 | Need n > 3*maxlag timepoints; lag 1-2 captures most direct regulation |
| dynGENIE3 threshold | Top 1000-5000 edges | 10-50 edges per TF for focused networks |
| DBN bootstrap R | 200 | Standard for moderate datasets; increase to 500 for noisy data |
| DBN strength cutoff | 0.7 | Edge in >70% of bootstraps; conservative but reliable |
| Granger FDR | q < 0.05 | Standard multiple testing threshold |

## Related Skills

- gene-regulatory-networks/coexpression-networks - Static co-expression networks
- gene-regulatory-networks/scenic-regulons - Single-cell regulon inference with pySCENIC
- gene-regulatory-networks/differential-networks - Condition-specific network comparison
- data-visualization/network-visualization - Network plotting with NetworkX and Cytoscape

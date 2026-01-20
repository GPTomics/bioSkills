# Similarity Network Fusion Usage Guide

## Overview

SNF integrates multiple data types by constructing patient similarity networks for each omics layer and fusing them into a single network that captures shared and complementary information.

## When to Use SNF

- **Patient stratification** - Identify disease subtypes
- **Prognosis prediction** - Survival-associated clusters
- **Heterogeneous data** - Different scales/distributions
- **Network-based analysis** - Sample relationships matter

## Key Parameters

### K (Neighbors)
Number of neighbors for affinity calculation. Typical: 10-30.

### alpha (Sigma scaling)
Controls kernel width. Typical: 0.3-0.8. Higher = broader similarity.

### t (Iterations)
Number of fusion iterations. Usually 10-20 sufficient.

## Algorithm Overview

1. **Distance matrices**: Compute pairwise sample distances per omics
2. **Affinity matrices**: Convert distances to similarities (Gaussian kernel)
3. **Network fusion**: Iteratively update networks using cross-view information
4. **Clustering**: Spectral clustering on fused network

## Data Preparation

- **Normalize** each omics separately (z-score, quantile)
- **Filter** low-variance features
- **Handle missing values** (impute or remove)
- **Match samples** across all omics

## Cluster Number Selection

- Use `estimateNumberOfClustersGivenGraph()`
- Consider biological interpretability
- Check stability across parameter choices

## Validation

- **NMI**: Compare to known labels if available
- **Survival**: Kaplan-Meier + log-rank test
- **Silhouette**: Cluster compactness

## References

- SNF paper: doi:10.1038/nmeth.2810
- SNFtool: https://cran.r-project.org/package=SNFtool

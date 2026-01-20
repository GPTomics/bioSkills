# Similarity Network Fusion - Usage Guide

## Overview
SNF integrates multiple data types by constructing patient similarity networks for each omics layer and fusing them into a single network that captures shared and complementary information.

## Prerequisites
```r
install.packages('SNFtool')
```
```python
pip install snfpy
```

## Quick Start
Tell your AI agent what you want to do:
- "Identify patient subtypes by fusing my RNA-seq and methylation data"
- "Build a similarity network from multiple omics to stratify patients"

## Example Prompts
### Patient Stratification
> "Use SNF to fuse my expression and CNV data to identify cancer subtypes"

> "Find patient clusters by integrating RNA, protein, and methylation similarity networks"

### Survival Analysis
> "Identify subtypes from multi-omics data and test for survival differences"

> "Stratify patients using SNF and generate Kaplan-Meier curves per cluster"

### Parameter Tuning
> "Run SNF with different K values to check cluster stability"

> "Estimate the optimal number of clusters from my fused similarity network"

## What the Agent Will Do
1. Load multi-omics data matrices with matched samples
2. Compute pairwise distance matrices per omics layer
3. Convert distances to affinity matrices using Gaussian kernel
4. Fuse networks using iterative SNF algorithm
5. Estimate optimal cluster number
6. Perform spectral clustering on fused network
7. Validate clusters (silhouette, survival, known labels)

## Key Parameters

| Parameter | Description | Typical Range |
|-----------|-------------|---------------|
| K | Number of neighbors for affinity | 10-30 |
| alpha | Kernel width scaling | 0.3-0.8 |
| t | Fusion iterations | 10-20 |

## Tips
- Normalize each omics separately (z-score or quantile) before computing distances
- Filter low-variance features to improve signal
- All samples must be present across all omics; impute or remove missing
- Use `estimateNumberOfClustersGivenGraph()` to choose cluster number
- Validate clusters with survival analysis if clinical data available
- Check cluster stability by varying K and alpha parameters

## References
- SNF paper: doi:10.1038/nmeth.2810
- SNFtool: https://cran.r-project.org/package=SNFtool

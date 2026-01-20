# mixOmics Analysis Usage Guide

## Overview

mixOmics provides multivariate methods for multi-omics integration, including both supervised (DIABLO, sPLS-DA) and unsupervised (sPLS, sPCA) approaches.

## Method Selection

| Method | Supervised | Blocks | Use Case |
|--------|------------|--------|----------|
| sPCA | No | 1 | Dimension reduction |
| sPLS | No | 2 | Find correlated features |
| sPLS-DA | Yes | 1 | Classify with feature selection |
| DIABLO | Yes | 2+ | Multi-omics classification |
| MINT | Yes | 1 | Multi-study integration |

## When to Use mixOmics

- **Classification** - Predict sample groups
- **Biomarker discovery** - Sparse feature selection
- **Multi-study** - Combine datasets with MINT
- **Correlation analysis** - Find co-varying features

## vs MOFA2

- **MOFA2**: Unsupervised, handles missing data, factor-based
- **mixOmics**: Supervised options, sparse selection, classification

## Key Parameters

### keepX
Number of features to select per component. Tune with cross-validation.

### ncomp
Number of latent components. Usually 2-5, check performance.

### design
For DIABLO: correlation structure between blocks (0 = independent, 1 = fully correlated).

## Data Preparation

1. Pre-filter low-variance features
2. Log-transform count data
3. Center and scale features
4. Remove outlier samples

## References

- mixOmics: doi:10.1371/journal.pcbi.1005752
- DIABLO: doi:10.1093/bioinformatics/bty1054
- Website: http://mixomics.org/

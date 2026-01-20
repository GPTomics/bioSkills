# MOFA2 Integration Usage Guide

## Overview

MOFA2 (Multi-Omics Factor Analysis v2) is an unsupervised method for integrating multiple omics layers. It decomposes the data into latent factors that explain shared and view-specific variation.

## When to Use MOFA2

- **Exploratory analysis** - Discover patterns across omics
- **Dimensionality reduction** - Multi-omics analog to PCA
- **Missing data** - MOFA handles partial overlap gracefully
- **No prior labels** - Unsupervised discovery

## Key Concepts

### Views
Different omics modalities (RNA, protein, methylation, etc.)

### Factors
Latent variables capturing sources of variation

### Weights
Feature loadings on each factor (for interpretation)

### Groups
Different experimental conditions or batches

## Data Preparation

1. **Normalization**: Pre-normalize each view separately
2. **Feature selection**: Select variable features per view (1000-5000)
3. **Centering**: Center features (mean = 0)
4. **Scaling**: Optional scaling (recommended for different units)

## Choosing Number of Factors

- Start with 15-25 factors
- MOFA automatically drops inactive factors
- Check variance explained to identify informative factors

## Interpretation Strategy

1. Examine variance explained per view per factor
2. Identify factors with shared vs view-specific variance
3. Extract top weighted features per factor
4. Run pathway enrichment on feature sets
5. Correlate factors with sample metadata

## References

- MOFA2 paper: doi:10.1186/s13059-020-02015-1
- MOFA+ (single-cell): doi:10.1038/s41592-019-0423-z
- Vignette: https://biofam.github.io/MOFA2/

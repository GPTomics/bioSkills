# Protein Quantification

## Overview

Protein quantification converts raw mass spectrometry signals into abundance estimates for statistical analysis.

## Quantification Strategies

### Label-Free Quantification (LFQ)
- **Intensity-based**: Sum of peptide intensities (MaxLFQ)
- **Spectral counting**: Number of spectra per protein
- **Pros**: No labeling required, unlimited samples
- **Cons**: More missing values, requires careful normalization

### Isobaric Labeling (TMT/iTRAQ)
- **TMT**: 6-plex, 10-plex, 11-plex, 16-plex, 18-plex
- **iTRAQ**: 4-plex, 8-plex
- **Pros**: Multiplex samples, reduced missing values
- **Cons**: Ratio compression, batch effects between plexes

### Metabolic Labeling (SILAC)
- Incorporate heavy amino acids during cell culture
- **Pros**: Most accurate ratios
- **Cons**: Limited to cell culture, max 3-plex

## Normalization Methods

| Method | Description | Use Case |
|--------|-------------|----------|
| Median centering | Shift to common median | General purpose |
| Quantile | Force identical distributions | Strong batch effects |
| LOESS | Local regression | Non-linear effects |
| VSN | Variance stabilization | Heteroscedastic data |

## Missing Value Handling

| Type | Method |
|------|--------|
| MCAR | Mean/median imputation |
| MAR | KNN imputation |
| MNAR (low abundance) | MinDet, MinProb, Perseus left-censored |

## Quality Metrics

- **CV per protein**: Coefficient of variation across replicates
- **Correlation**: Pearson/Spearman between replicates
- **PCA**: Check for batch effects and outliers
- **Missing value pattern**: Should be random, not systematic

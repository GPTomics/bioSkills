# Multi-omics Integration Pipeline Usage Guide

## Overview

This workflow integrates multiple molecular data types (transcriptomics, proteomics, metabolomics, etc.) to discover shared biological signals and cross-modal biomarkers.

## Prerequisites

```r
BiocManager::install(c('MOFA2', 'mixOmics'))
install.packages(c('SNFtool', 'pheatmap'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Integrate my transcriptomics and proteomics data"
- "Run MOFA2 on my multi-omics dataset"
- "Find shared factors across my omics modalities"

## Example Prompts

### Basic Integration
> "I have RNA-seq and proteomics from the same samples, integrate them with MOFA2"

> "Find shared biological signals across my multi-omics data"

### Supervised Analysis
> "Use DIABLO to find multi-omics biomarkers that predict treatment response"

> "Run supervised multi-omics analysis with patient outcomes"

### Patient Stratification
> "Cluster patients using multi-omics data with SNF"

> "Find patient subtypes from my integrated omics data"

## When to Use This Pipeline

- Studies with multiple omics measurements
- Multi-platform biomarker discovery
- Patient stratification across modalities
- Understanding disease mechanisms
- Drug response prediction

## Required Inputs

1. **Normalized data matrices** - One per modality (samples × features)
2. **Sample metadata** - Conditions, outcomes, covariates
3. **Feature annotations** - Gene names, protein IDs, metabolite names

## Data Requirements

### Samples
- Common samples across modalities (or at least significant overlap)
- Minimum 20-30 samples for MOFA, more for mixOmics
- Biological replicates preferred

### Features
- Pre-normalized data (log2, z-score, etc.)
- Similar scales across modalities
- Feature IDs that can be mapped (genes, proteins, metabolites)

## Integration Methods

### MOFA2 (Recommended)
**Best for:** Unsupervised exploration, factor discovery

- Handles missing views/samples
- Identifies shared vs. unique variation
- Easy interpretation via factor weights
- GPU acceleration available

### mixOmics DIABLO
**Best for:** Supervised classification, biomarker selection

- Requires outcome variable
- Feature selection built-in
- Cross-validation for variable tuning
- Sparse solutions

### SNF (Similarity Network Fusion)
**Best for:** Patient clustering, network-based integration

- Creates fused similarity network
- Good for patient stratification
- No feature selection
- Works with any data type

## Pipeline Steps

### 1. Data Preprocessing
- Normalize each modality separately
- Log transform if needed
- Handle missing values

### 2. Sample Harmonization
- Match samples across modalities
- Verify sample IDs are consistent
- Remove samples missing from all modalities

### 3. Feature Selection
- Select most variable features
- Optional: pre-filter by relevance
- Balance features across modalities

### 4. Integration
- Run MOFA/DIABLO/SNF
- Tune parameters (components, neighbors)
- Assess convergence

### 5. Interpretation
- Analyze factors/components
- Identify top features per factor
- Pathway enrichment on feature sets

## Parameter Guidelines

### MOFA2
| Parameter | Recommended | Notes |
|-----------|-------------|-------|
| num_factors | 10-15 | Start high, prune later |
| scale_views | TRUE | Equalizes modality contributions |
| maxiter | 1000 | Increase if not converged |

### mixOmics DIABLO
| Parameter | Recommended | Notes |
|-----------|-------------|-------|
| ncomp | 2-5 | Tune via CV |
| keepX | 10-100 | Tune per modality |
| design | 0.1 | Low correlation between blocks |

### SNF
| Parameter | Recommended | Notes |
|-----------|-------------|-------|
| K (neighbors) | 10-30 | ~10-20% of samples |
| t (iterations) | 20 | Usually sufficient |

## Common Issues

### Few common samples
- Check sample ID formatting
- Consider imputation for missing views
- MOFA handles missing better than DIABLO

### Factors dominated by one view
- Increase scale_views strength
- Balance feature numbers
- Check for outliers in dominant view

### No meaningful factors
- Check data quality
- Try different feature selection
- Verify biological signal exists

### DIABLO overfitting
- Reduce number of features
- Use more CV folds
- Check sample size vs features

## Output Files

| File | Description |
|------|-------------|
| mofa_model.hdf5 | Trained MOFA model |
| mofa_factor_values.csv | Sample factor scores |
| mofa_weights.csv | Feature weights per factor |
| variance_explained.png | Factor variance decomposition |
| factor_scatter.png | Sample visualization |
| top_weights_*.png | Important features |

## Interpretation Guide

### Factor Variance Explained
- Factors explaining >5% variance are meaningful
- Shared factors: affect multiple modalities
- Unique factors: one modality only

### Feature Weights
- High absolute weight = important feature
- Sign indicates direction of association
- Compare top features across modalities

### Biological Interpretation
1. Extract top genes per factor
2. Run pathway enrichment
3. Connect to clinical variables
4. Validate key features

## Tips

- **Sample overlap**: Ensure samples are matched across modalities (same patients/conditions)
- **Normalization**: Pre-normalize each modality separately before integration
- **Missing views**: MOFA2 handles missing views better than DIABLO
- **Feature selection**: Select top variable features per modality to reduce noise
- **Interpretation**: Run pathway enrichment on top-weighted features per factor

## References

- MOFA2: doi:10.1186/s13059-020-02015-1
- mixOmics: doi:10.1371/journal.pcbi.1005752
- SNF: doi:10.1038/nmeth.2810

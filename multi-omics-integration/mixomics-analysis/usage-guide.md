# mixOmics Analysis - Usage Guide

## Overview
mixOmics provides multivariate methods for multi-omics integration, including both supervised (DIABLO, sPLS-DA) and unsupervised (sPLS, sPCA) approaches with sparse feature selection.

## Prerequisites
```r
BiocManager::install('mixOmics')
```

## Quick Start
Tell your AI agent what you want to do:
- "Classify my samples using RNA and protein data with DIABLO"
- "Find correlated features between transcriptomics and metabolomics"

## Example Prompts
### Classification
> "Use DIABLO to classify tumor vs normal samples using my RNA-seq and proteomics data"

> "Run sPLS-DA on my expression matrix to identify biomarkers distinguishing treatment groups"

### Feature Selection
> "Find the most discriminative features across my three omics layers for predicting response"

> "Select sparse biomarker panels from each omics view using DIABLO"

### Correlation Analysis
> "Identify genes and metabolites that co-vary across my samples using sPLS"

> "Find correlated protein-metabolite pairs driving the separation between groups"

### Multi-Study
> "Combine RNA-seq from three studies using MINT to find robust biomarkers"

## What the Agent Will Do
1. Load and preprocess multi-omics matrices
2. Select appropriate method (sPLS, sPLS-DA, DIABLO, MINT)
3. Tune keepX parameters using cross-validation
4. Build the integration model
5. Extract selected features per component
6. Generate visualization plots (sample plots, loading plots, circos)
7. Assess model performance if supervised

## Method Selection Guide

| Method | Supervised | Blocks | Use Case |
|--------|------------|--------|----------|
| sPCA | No | 1 | Dimension reduction |
| sPLS | No | 2 | Find correlated features |
| sPLS-DA | Yes | 1 | Classify with feature selection |
| DIABLO | Yes | 2+ | Multi-omics classification |
| MINT | Yes | 1 | Multi-study integration |

## Tips
- Use MOFA2 for unsupervised discovery; mixOmics excels at supervised classification
- Tune keepX (features per component) with cross-validation for best performance
- ncomp usually 2-5; check classification performance to choose
- For DIABLO, design matrix controls correlation structure between blocks (0=independent, 1=correlated)
- Pre-filter low-variance features and log-transform count data before analysis

## References
- mixOmics: doi:10.1371/journal.pcbi.1005752
- DIABLO: doi:10.1093/bioinformatics/bty1054
- Website: http://mixomics.org/

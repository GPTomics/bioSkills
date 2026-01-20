# MOFA2 Integration - Usage Guide

## Overview
MOFA2 (Multi-Omics Factor Analysis v2) is an unsupervised method for integrating multiple omics layers. It decomposes the data into latent factors that explain shared and view-specific variation.

## Prerequisites
```bash
pip install mofapy2
```
```r
BiocManager::install('MOFA2')
```

## Quick Start
Tell your AI agent what you want to do:
- "Integrate my RNA-seq and proteomics data to find shared patterns"
- "Run MOFA2 on my multi-omics experiment to identify latent factors"

## Example Prompts
### Basic Integration
> "Run MOFA2 on my expression and methylation matrices to discover shared factors"

> "Create a MOFA model from my scRNA-seq and scATAC-seq data"

### Interpretation
> "Which factors explain the most variance across my omics views?"

> "Extract the top features driving factor 1 and run pathway enrichment"

### Visualization
> "Plot the variance explained per factor across RNA and protein views"

> "Correlate MOFA factors with my sample metadata to find clinical associations"

## What the Agent Will Do
1. Load and validate your multi-omics data matrices
2. Create MOFA object with appropriate views and groups
3. Set model options (number of factors, convergence criteria)
4. Train the MOFA model
5. Analyze variance explained per factor per view
6. Extract top-weighted features for biological interpretation
7. Generate visualization plots

## Key Concepts

### Views
Different omics modalities (RNA, protein, methylation, etc.)

### Factors
Latent variables capturing sources of variation

### Weights
Feature loadings on each factor (for interpretation)

### Groups
Different experimental conditions or batches

## Tips
- Start with 15-25 factors; MOFA automatically drops inactive factors
- Pre-normalize each view separately before integration
- Select 1000-5000 variable features per view for best results
- Center features (mean = 0) and optionally scale when units differ
- Check variance explained to identify informative factors
- MOFA handles partial sample overlap gracefully (missing data OK)

## References
- MOFA2 paper: doi:10.1186/s13059-020-02015-1
- MOFA+ (single-cell): doi:10.1038/s41592-019-0423-z
- Vignette: https://biofam.github.io/MOFA2/

# Proteomics QC - Usage Guide

## Overview
Assess data quality before statistical analysis through sample correlation, missing value patterns, batch effects, and replicate reproducibility.

## Prerequisites
```bash
pip install numpy pandas scipy scikit-learn matplotlib seaborn
# R packages: BiocManager::install(c("limma", "sva", "pcaMethods"))
```

## Quick Start
Tell your AI agent what you want to do:
- "Run QC on my protein matrix and identify outlier samples"
- "Check for batch effects in my proteomics data"
- "Generate a QC report with correlation heatmap and PCA"

## Example Prompts

### Sample Quality
> "Calculate the correlation matrix between samples and identify any with correlation < 0.9"

> "Count proteins identified per sample and flag samples with significantly fewer"

> "Check median intensity per sample for outliers"

### Missing Values
> "Analyze the missing value pattern - is it random or systematic?"

> "Create a heatmap showing missing values by sample and protein"

> "Calculate the percentage of missing values per sample and per protein"

### Batch Effects
> "Run PCA and color by batch to check for batch effects"

> "Test if PC1 or PC2 separates by run date rather than biology"

> "Recommend whether ComBat or SVA is needed for my data"

### Replicate Reproducibility
> "Calculate CV per protein across technical replicates"

> "Check if biological replicates cluster together in PCA"

> "Compare intensity distributions between replicate groups"

### QC Report
> "Generate a comprehensive QC report with all standard metrics"

> "Create QC plots: correlation heatmap, PCA, intensity boxplots, missing values"

## What the Agent Will Do
1. Load and log2-transform protein matrix
2. Calculate sample-level metrics (proteins, missing %, median intensity)
3. Compute correlation matrix and flag outliers
4. Run PCA and check for batch effects
5. Analyze missing value patterns
6. Calculate replicate CV
7. Generate QC report and visualizations

## Key QC Metrics

### Sample-Level
| Metric | Good | Warning | Action |
|--------|------|---------|--------|
| Proteins identified | >2000 | <1500 | Check MS performance |
| Missing values | <30% | >50% | Consider exclusion |
| Replicate correlation | >0.95 | <0.90 | Check sample prep |
| Median CV | <20% tech | >40% | Check variability source |

### Batch Effect Signs
- PC1/PC2 separates by batch, not biology
- Samples cluster by run date
- Different intensity distributions by batch

## Tips
- Always run QC before statistical analysis
- Document and justify any sample exclusions
- Include batch as covariate if effects detected
- Consider sensitivity analysis with/without borderline samples
- Technical CV <20%, biological CV <40% are typical thresholds

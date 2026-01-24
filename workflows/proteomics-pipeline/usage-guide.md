# Proteomics Pipeline Usage Guide

## Overview

End-to-end workflow for label-free proteomics analysis from MaxQuant/DIA-NN output to differential protein abundance.

## Prerequisites

```r
BiocManager::install(c('limma', 'DEP', 'MSstats'))
install.packages(c('pheatmap', 'ggplot2'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Run the proteomics pipeline on my MaxQuant output"
- "Find differentially expressed proteins between conditions"
- "Process my DIA-NN results and run differential analysis"

## Example Prompts

### Basic Analysis
> "I have proteinGroups.txt from MaxQuant, run the full pipeline"

> "Normalize my proteomics data and find differential proteins"

### QC and Preprocessing
> "Check sample quality with PCA and correlation heatmap"

> "Impute missing values in my proteomics data using MinProb"

### Differential Analysis
> "Run limma to find proteins changed between treatment and control"

> "Use MSstats for differential analysis with my peptide-level data"

## Pipeline Stages

### 1. Data Import
- Load proteinGroups.txt (MaxQuant) or report.tsv (DIA-NN)
- Filter contaminants and decoys
- Extract intensity columns

### 2. Transformation
- Replace 0 with NA
- Log2 transform
- Median centering normalization

### 3. Filtering
- Remove proteins with >50% missing values
- Filter low-variance proteins (optional)

### 4. Imputation
- MinProb: Left-censored Gaussian (for MNAR)
- KNN: K-nearest neighbors (for MAR)
- Perseus-style: Downshifted Gaussian

### 5. Quality Control
- PCA: Check replicate clustering
- Correlation heatmap: Sample similarity
- Missing value patterns: Random or systematic

### 6. Differential Analysis
- limma: Empirical Bayes moderated t-test
- MSstats: Mixed-effects models
- DEP: Complete workflow package

### 7. Output
- Differential proteins table
- Volcano plot
- Heatmap of significant proteins

## Input Requirements

### MaxQuant Output
```
proteinGroups.txt  # Protein-level quantification
evidence.txt       # Peptide-level (for MSstats)
annotation.csv     # Sample metadata
```

### Sample Annotation
```csv
sample,condition,replicate
Sample1,Control,1
Sample2,Control,2
Sample3,Treatment,1
Sample4,Treatment,2
```

## Expected Outputs

| File | Description |
|------|-------------|
| differential_proteins.csv | All proteins with statistics |
| volcano_plot.pdf | Log2FC vs -log10(p-value) |
| pca_plot.pdf | Sample clustering |
| heatmap.pdf | Significant proteins |

## Typical Results

- 2000-5000 quantified proteins (cell lysate)
- 50-500 differential proteins (10%)
- Fold changes typically 1.5-4x

## Tips

- **Missing values**: Use MinProb for MNAR (abundance-dependent), KNN for MAR
- **Normalization**: Median centering is standard; quantile if distributions vary
- **Filtering**: Remove proteins with >50% missing values before imputation
- **Replicates**: Minimum 3 biological replicates per condition
- **Contaminants**: Always filter MaxQuant contaminants and reverse sequences

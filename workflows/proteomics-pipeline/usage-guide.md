# Proteomics Pipeline

## Overview

End-to-end workflow for label-free proteomics analysis from MaxQuant/DIA-NN output to differential protein abundance.

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

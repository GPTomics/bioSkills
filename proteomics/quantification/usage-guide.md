# Quantification - Usage Guide

## Overview
Convert raw mass spectrometry signals into protein abundance estimates using label-free, isobaric, or metabolic labeling strategies.

## Prerequisites
```bash
pip install numpy pandas scipy
# R packages: BiocManager::install(c("MSstats", "DEP", "MSnbase"))
```

## Quick Start
Tell your AI agent what you want to do:
- "Normalize my MaxLFQ intensities using median centering"
- "Process TMT reporter ion intensities from my experiment"
- "Impute missing values using KNN for my protein matrix"

## Example Prompts

### Normalization
> "Apply median centering normalization to my protein intensity matrix"

> "Use quantile normalization to correct for batch effects between runs"

> "Normalize my TMT data using the internal reference channel"

### Missing Value Handling
> "Impute missing values using KNN for MAR pattern and MinProb for MNAR"

> "Filter proteins with more than 50% missing values, then impute the rest"

> "Analyze the missing value pattern and recommend an imputation strategy"

### Label-Free Quantification
> "Calculate MaxLFQ intensities from peptide-level data"

> "Summarize peptide intensities to protein level using top3 method"

### TMT/iTRAQ Processing
> "Extract TMT reporter ion intensities and correct for isotope impurity"

> "Normalize across TMT plexes using the bridge channel"

## What the Agent Will Do
1. Load protein/peptide intensity matrix
2. Log2-transform raw intensities
3. Apply appropriate normalization method
4. Identify missing value pattern (MCAR/MAR/MNAR)
5. Impute missing values with suitable method
6. Generate QC metrics (CV, correlation)

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
| MNAR (low abundance) | MinDet, MinProb, left-censored |

## Tips
- Always log2-transform before normalization
- Check CV across replicates (technical <20%, biological <40%)
- Use PCA to verify normalization removed batch effects
- Document imputation method - it affects downstream statistics

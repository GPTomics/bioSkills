# Normalization and QC - Usage Guide

## Overview

Proper normalization and QC are critical for metabolomics data quality. Batch effects and instrumental drift can dominate biological signals if not corrected.

## Prerequisites

```bash
# R packages
install.packages(c("statTarget", "MetNorm"))
BiocManager::install("sva")  # ComBat

# Python
pip install sklearn pandas numpy
```

## Quick Start

Tell your AI agent what you want to do:
- "Apply QC-RSC normalization to correct for instrument drift"
- "Check QC sample clustering and calculate feature RSDs"

## Example Prompts

### Quality Assessment
> "Calculate RSD for each feature across QC samples and flag features with RSD > 30%"
> "Run PCA to check if QC samples cluster together"
> "Check for batch effects in my metabolomics data"

### Normalization
> "Apply QC-RSC drift correction using my pooled QC samples"
> "Use PQN normalization across all samples"
> "Run ComBat batch correction with sample groups as covariates"

### Transformation
> "Log2 transform the intensity data after normalization"
> "Apply Pareto scaling before multivariate analysis"

### Missing Values
> "Filter features with more than 30% missing values"
> "Impute missing values using KNN within each sample group"

## What the Agent Will Do

1. Assess data quality (RSDs, missing values, TICs)
2. Check for batch effects via PCA
3. Apply drift correction using QC samples
4. Normalize sample loading differences
5. Transform and scale for analysis
6. Export normalized feature table

## Tips

- Include pooled QC samples every 10 biological samples
- QC RSD < 30% is acceptable (< 20% for targeted assays)
- QC samples should cluster in PCA; separation indicates problems
- Always log transform MS intensity data before statistics
- Pareto scaling reduces influence of high-abundance features

## Normalization Methods

| Method | Use When |
|--------|----------|
| QC-RSC | QC samples available, instrumental drift |
| TIC | No QC, simple correction |
| PQN | NMR or general purpose |
| ComBat | Known batch structure |

## Quality Thresholds

| Metric | Acceptable | Concerning |
|--------|------------|------------|
| QC RSD | < 30% | > 50% |
| Missing values | < 20% | > 50% |
| QC clustering | Tight | Dispersed |

## References

- QC-RSC: doi:10.1007/s11306-016-1030-9
- PQN: doi:10.1021/ac051632c
- ComBat: doi:10.1093/biostatistics/kxj037

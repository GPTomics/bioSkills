# Normalization and QC Usage Guide

## Overview

Proper normalization and QC are critical for metabolomics data quality. Batch effects and instrumental drift can dominate biological signals.

## Normalization Methods

| Method | Use When |
|--------|----------|
| QC-RSC | QC samples available, instrumental drift |
| TIC | No QC, simple correction |
| PQN | NMR or general purpose |
| ComBat | Known batch structure |

## QC Strategy

### Pooled QC Samples
- Mix of all samples
- Inject every 10 biological samples
- Use for drift correction and quality monitoring

### Quality Metrics
- QC RSD <30% (20% for targeted)
- QC samples cluster in PCA
- Coefficient of variation across QCs

## Transformation Guidelines

### Log Transformation
- Always for MS intensity data
- Stabilizes variance
- log2 or log10

### Scaling
- **Pareto**: Reduces influence of high-abundance features
- **Auto (z-score)**: Equal weight to all features
- **Range**: Scale to 0-1

## Missing Value Strategy

1. Filter features with >20-50% missing
2. Impute by group (biological replicates)
3. Methods: KNN, minimum value, BPCA

## Red Flags

- QC RSD >50% (instrument issues)
- QC samples don't cluster
- Strong batch effect in PCA
- >50% missing values

## References

- QC-RSC: doi:10.1007/s11306-016-1030-9
- PQN: doi:10.1021/ac051632c

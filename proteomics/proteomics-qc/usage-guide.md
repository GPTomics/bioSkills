# Proteomics Quality Control

## Overview

Quality control ensures data quality before statistical analysis. Poor QC can lead to false positives/negatives and irreproducible results.

## Key QC Metrics

### Sample-Level
| Metric | Good | Warning | Action |
|--------|------|---------|--------|
| Proteins identified | >2000 | <1500 | Check MS performance |
| Missing values | <30% | >50% | Consider exclusion |
| Replicate correlation | >0.95 | <0.90 | Check sample prep |
| Median CV | <20% tech, <40% bio | Higher | Check variability source |

### Protein-Level
| Metric | Threshold | Action |
|--------|-----------|--------|
| Missing per protein | <50% | Filter proteins |
| Intensity distribution | Normal (log2) | Check normalization |

## QC Workflow

1. **Load and filter** data (contaminants, reverse)
2. **Log2 transform** intensities
3. **Check distributions** - should be approximately normal
4. **Correlation analysis** - replicates should cluster
5. **PCA** - check for batch effects
6. **Missing value analysis** - identify patterns
7. **CV calculation** - assess technical variation
8. **Generate report** - document QC metrics

## Batch Effects

Signs of batch effects:
- PC1/PC2 separates by batch, not biology
- Samples cluster by run date
- Different intensity distributions by batch

Remediation:
- ComBat for known batches
- SVA for unknown batches
- Include batch as covariate

## Missing Value Types

| Type | Pattern | Cause | Handling |
|------|---------|-------|----------|
| MCAR | Random | Technical | Any imputation |
| MAR | Covariate-dependent | Sample-related | KNN, regression |
| MNAR | Low abundance | Below detection | MinProb, left-censored |

## Outlier Detection

Methods:
- PCA: samples far from cluster
- Correlation: low correlation with replicates
- Protein count: significantly fewer proteins
- Median intensity: very different from others

Action:
- Document reason for exclusion
- Consider sensitivity analysis with/without

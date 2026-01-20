# Data Harmonization Usage Guide

## Overview

Before multi-omics integration, data must be harmonized to ensure compatibility across data types. This includes normalization, batch correction, feature alignment, and handling missing values.

## Preprocessing Steps

### 1. Quality Control (per assay)
- Filter low-quality samples
- Remove outliers
- Check for batch effects

### 2. Normalization (assay-specific)
| Data Type | Method |
|-----------|--------|
| RNA-seq counts | VST, rlog, TMM |
| Proteomics intensity | Log2 + median centering |
| Methylation beta | M-value transform |
| Metabolomics | Log + pareto scaling |

### 3. Batch Correction
- ComBat (parametric or non-parametric)
- limma removeBatchEffect
- SVA for unknown batches

### 4. Feature Alignment
- Map features to common namespace (gene symbols)
- Aggregate if needed (proteins to genes)

### 5. Missing Values
- Filter features with >30-50% missing
- Impute remaining (MinProb, KNN, etc.)

### 6. Scaling
- Z-score (mean=0, sd=1)
- Or quantile normalization across assays

## Common Issues

### Different Sample Sets
Not all samples have all omics. Options:
- Use only complete samples
- Impute at sample level
- Use methods tolerant to missing views (MOFA2)

### Different Feature Spaces
RNA has genes, proteomics has proteins. Options:
- Map to common genes
- Keep separate (integration handles it)
- Feature-level integration only on matched

### Scale Differences
Expression values vs methylation betas. Always:
- Normalize within assay first
- Scale/center before integration

## References

- MultiAssayExperiment: doi:10.1158/0008-5472.CAN-17-0344
- ComBat: doi:10.1093/biostatistics/kxj037

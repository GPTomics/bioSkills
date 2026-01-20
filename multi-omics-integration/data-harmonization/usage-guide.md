# Data Harmonization - Usage Guide

## Overview
Before multi-omics integration, data must be harmonized to ensure compatibility across data types. This includes normalization, batch correction, feature alignment, and handling missing values.

## Prerequisites
```r
BiocManager::install(c('MultiAssayExperiment', 'sva', 'limma'))
```
```python
pip install combat numpy pandas
```

## Quick Start
Tell your AI agent what you want to do:
- "Prepare my RNA-seq and proteomics data for MOFA2 integration"
- "Batch correct and normalize my multi-omics experiment"

## Example Prompts
### Normalization
> "Normalize my RNA-seq counts with VST and proteomics intensities with log2 median centering"

> "Apply appropriate normalization to each of my omics layers before integration"

### Batch Correction
> "Remove batch effects from my expression data using ComBat before integrating with proteomics"

> "Apply limma removeBatchEffect to correct for sequencing batch across samples"

### Sample Matching
> "Match samples across my RNA, protein, and methylation datasets by sample ID"

> "Identify which samples have complete data across all my omics assays"

### Feature Alignment
> "Map protein IDs to gene symbols to align with my RNA-seq data"

> "Aggregate protein-level measurements to gene-level for integration"

### Missing Data
> "Impute missing values in my proteomics data using MinProb before integration"

> "Filter features with more than 30% missing values across samples"

## What the Agent Will Do
1. Assess data quality per omics layer
2. Apply assay-specific normalization
3. Detect and correct batch effects
4. Align feature identifiers across omics
5. Handle missing values (filter or impute)
6. Scale features for integration
7. Create MultiAssayExperiment object if using R

## Normalization by Data Type

| Data Type | Recommended Method |
|-----------|-------------------|
| RNA-seq counts | VST, rlog, or TMM |
| Proteomics intensity | Log2 + median centering |
| Methylation beta | M-value transform |
| Metabolomics | Log + pareto scaling |

## Tips
- Always normalize within assay before scaling across assays
- Use ComBat for known batches; SVA for unknown confounders
- Map features to common namespace (gene symbols) when possible
- Filter features with >30-50% missing before imputation
- Methods like MOFA2 tolerate missing views; SNF requires complete overlap
- Z-score scaling (mean=0, sd=1) is standard before most integration methods

## References
- MultiAssayExperiment: doi:10.1158/0008-5472.CAN-17-0344
- ComBat: doi:10.1093/biostatistics/kxj037

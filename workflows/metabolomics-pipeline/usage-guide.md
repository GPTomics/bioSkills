# Metabolomics Pipeline Usage Guide

## Overview

This workflow processes raw mass spectrometry data through peak detection, alignment, normalization, statistical analysis, and pathway interpretation.

## Prerequisites

```r
BiocManager::install(c('xcms', 'CAMERA', 'MetaboAnalystR'))
install.packages(c('metablastr', 'pheatmap'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Run the metabolomics pipeline on my mzML files"
- "Process my LC-MS data and find differential metabolites"
- "Analyze my lipidomics experiment"

## Example Prompts

### Basic Analysis
> "I have mzML files from an untargeted metabolomics study, run the full pipeline"

> "Process my LC-MS/MS data with XCMS and run differential analysis"

### Normalization and QC
> "Apply QC-based batch correction to my metabolomics data"

> "Normalize my metabolomics data and check sample quality with PCA"

### Pathway Analysis
> "Find enriched metabolic pathways in my differential metabolites"

> "Annotate my significant features against HMDB and run pathway enrichment"

## When to Use This Pipeline

- Untargeted metabolomics studies
- LC-MS/MS metabolite profiling
- Lipidomics analysis
- Metabolic biomarker discovery
- Treatment response studies

## Required Inputs

1. **Raw MS data** - mzML or mzXML format (converted from vendor formats)
2. **Sample metadata** - CSV with sample names, conditions, batches
3. **QC samples** - Pooled QC samples recommended

## Sample Metadata Format

```csv
sample,condition,batch,injection_order
Sample1.mzML,Control,1,1
Sample2.mzML,Control,1,2
QC1.mzML,QC,1,3
Sample3.mzML,Treatment,1,4
```

## Pipeline Steps

### 1. Peak Detection
- Identifies chromatographic peaks in each sample
- CentWave algorithm for LC-MS data
- Adjust peakwidth based on chromatography

### 2. Retention Time Alignment
- Corrects RT drift between samples
- Obiwarp or peak groups methods
- Essential for feature matching

### 3. Feature Grouping
- Groups peaks across samples into features
- Based on m/z and aligned RT
- minFraction controls stringency

### 4. Gap Filling
- Recovers missing values
- Integrates signal at expected locations
- Reduces false missing values

### 5. Normalization
- Corrects systematic variation
- Options: median, quantile, LOESS
- QC-based correction for batches

### 6. Statistical Analysis
- limma for differential analysis
- Handles missing values
- Multiple testing correction

### 7. Annotation
- Match m/z to databases (HMDB, KEGG, LipidMaps)
- Consider adducts and isotopes
- MS/MS matching for confidence

### 8. Pathway Analysis
- Map to KEGG pathways
- Over-representation analysis
- Metabolite set enrichment

## Parameter Guidelines

### Peak Detection (CentWave)
| Parameter | UPLC | Standard LC | GC-MS |
|-----------|------|-------------|-------|
| peakwidth | 5-30 | 10-60 | 2-10 |
| ppm | 15-25 | 25-50 | 10-20 |
| snthresh | 10 | 10 | 5 |

### Feature Grouping
| Parameter | Typical | Stringent |
|-----------|---------|-----------|
| bw | 5-10 | 2-3 |
| minFraction | 0.5 | 0.8 |
| binSize | 0.025 | 0.01 |

## Quality Control

### QC Sample Strategy
- Pool equal volumes from all samples
- Inject QC every 5-10 samples
- Use for batch correction and quality assessment

### QC Metrics
| Metric | Good | Acceptable | Poor |
|--------|------|------------|------|
| Features detected | >5000 | 2000-5000 | <2000 |
| QC CV | <20% | 20-30% | >30% |
| Blank ratio | >10x | 5-10x | <5x |

## Common Issues

### Few features detected
- Adjust peak detection parameters
- Check raw data quality
- Lower snthresh carefully

### Poor alignment
- Check for RT drift pattern
- Use more reference peaks
- Consider subset alignment

### High missing values
- Reduce minFraction
- Improve gap filling
- Check sample quality

### No significant features
- Check experimental design
- Consider effect sizes
- Adjust FDR threshold

## Output Files

| File | Description |
|------|-------------|
| normalized_feature_matrix.csv | Processed feature intensities |
| differential_metabolites.csv | Statistical results |
| qc_pca.png | PCA quality check |
| volcano_metabolites.png | Differential analysis plot |
| pathway_overview.png | Enriched pathways |

## Tips

- **QC samples**: Inject pooled QC every 5-10 samples for batch correction
- **Peak detection**: Adjust peakwidth based on your chromatography (UPLC: 5-30, standard LC: 10-60)
- **Missing values**: High missing values may indicate poor sample quality
- **Annotation confidence**: MS/MS matching provides higher confidence than m/z alone
- **mzML conversion**: Convert vendor files using ProteoWizard msConvert

## References

- XCMS: doi:10.1021/ac051437y
- MetaboAnalystR: doi:10.1093/bioinformatics/btaa123
- xcms3 workflow: doi:10.3390/metabo10120504

# methylKit Analysis - Usage Guide

## Overview
methylKit is an R/Bioconductor package for DNA methylation analysis, handling data import, quality control, normalization, and differential methylation analysis from bisulfite sequencing data.

## Prerequisites
```r
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')

BiocManager::install('methylKit')
```

## Quick Start
Tell your AI agent what you want to do:
- "Load my Bismark coverage files into methylKit and run QC"
- "Find differentially methylated CpGs between treatment and control"
- "Generate PCA and correlation plots for my methylation samples"

## Example Prompts
### Data Import
> "Read my Bismark coverage files into methylKit with sample metadata"

> "Import methylation data for 4 controls and 4 treated samples"

### Quality Control
> "Generate methylation and coverage statistics plots for all samples"

> "Show me PCA and sample correlation for my methylation data"

### Differential Analysis
> "Find differentially methylated CpGs with at least 25% difference and q < 0.01"

> "Run differential methylation analysis between tumor and normal samples"

> "Identify hypermethylated and hypomethylated CpGs separately"

### Filtering
> "Filter CpGs by coverage (minimum 10x) and normalize samples"

> "Unite samples requiring CpG coverage in at least 3 samples per group"

## What the Agent Will Do
1. Create sample metadata with file paths, sample IDs, and treatment groups
2. Import Bismark coverage files with methRead()
3. Generate QC plots (coverage stats, methylation stats)
4. Filter by coverage and normalize between samples
5. Unite samples to get common CpGs
6. Visualize sample relationships (PCA, correlation)
7. Run differential methylation analysis
8. Export significant differentially methylated CpGs

## Output Interpretation

| Column | Description |
|--------|-------------|
| chr, start, end | Genomic position |
| meth.diff | Methylation difference (%) |
| pvalue | Raw p-value |
| qvalue | FDR-adjusted p-value |

Positive meth.diff = hypermethylated in treatment
Negative meth.diff = hypomethylated in treatment

## Tips
- Use pipeline = 'bismarkCoverage' when reading Bismark .cov files
- Set destrand = TRUE in unite() to combine CpGs on both strands
- Typical filters: lo.count = 10 (minimum coverage), hi.perc = 99.9 (remove PCR artifacts)
- For memory issues, use save.db = TRUE for database-backed objects
- Use min.per.group in unite() if samples have variable coverage
- Overdispersion = 'MN' (multiplicative) is recommended for calculateDiffMeth()
- Common thresholds: difference = 25%, qvalue = 0.01

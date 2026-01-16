# Differential Binding Analysis with DiffBind - Usage Guide

## Overview

DiffBind identifies differentially bound regions between ChIP-seq conditions using read counts in consensus peaks. It applies normalization and statistical testing similar to RNA-seq differential expression analysis.

## When to Use This Skill

- You have ChIP-seq data from multiple conditions with replicates
- You want to find peaks that change between conditions
- You need statistical significance for differential binding

## Installation

```r
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')

BiocManager::install('DiffBind')
```

## Requirements

- Aligned BAM files for each sample
- Peak files from each sample (MACS2 narrowPeak/broadPeak)
- At least 2 replicates per condition (recommended: 3+)

## Sample Sheet Format

Create a CSV file with columns:
- SampleID: Unique sample name
- Condition: Control vs Treatment
- Replicate: 1, 2, 3...
- bamReads: Path to ChIP BAM
- Peaks: Path to peak file
- PeakCaller: macs (for narrowPeak)

## Basic Workflow

```r
library(DiffBind)

# 1. Load samples
dba_obj <- dba(sampleSheet = 'samples.csv')

# 2. Count reads
dba_obj <- dba.count(dba_obj, summits = 250)

# 3. Normalize
dba_obj <- dba.normalize(dba_obj)

# 4. Set contrast
dba_obj <- dba.contrast(dba_obj, categories = DBA_CONDITION)

# 5. Analyze
dba_obj <- dba.analyze(dba_obj)

# 6. Get results
results <- dba.report(dba_obj, th = 0.05)
```

## Common Issues

### Few Differential Peaks

- Ensure sufficient replicates (n >= 2)
- Check peak quality and overlap between replicates
- Try relaxing FDR threshold

### Memory Issues

- Use `summits = FALSE` to keep original peak widths
- Reduce number of peaks with `minOverlap`

## Output Interpretation

| Column | Description |
|--------|-------------|
| Conc | Mean concentration (log2) |
| Conc_control | Control mean |
| Conc_treatment | Treatment mean |
| Fold | Log2 fold change |
| p-value | Raw p-value |
| FDR | Adjusted p-value |

Positive Fold = gained in treatment
Negative Fold = lost in treatment

## Resources

- [DiffBind Bioconductor](https://bioconductor.org/packages/DiffBind/)
- [DiffBind Vignette](https://bioconductor.org/packages/release/bioc/vignettes/DiffBind/inst/doc/DiffBind.html)

# methylKit Analysis - Usage Guide

## Overview

methylKit is an R/Bioconductor package for analysis of DNA methylation data from bisulfite sequencing. It handles data import, quality control, normalization, and differential methylation analysis.

## When to Use This Skill

- You have Bismark coverage files from multiple samples
- You want to identify differentially methylated CpGs
- You need to visualize methylation patterns
- You want to compare methylation between conditions

## Installation

```r
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')

BiocManager::install('methylKit')
```

## Basic Workflow

### 1. Prepare Sample Information

Create a sample sheet:
```r
samples <- data.frame(
    file = c('ctrl1.cov.gz', 'ctrl2.cov.gz', 'treat1.cov.gz', 'treat2.cov.gz'),
    sample_id = c('ctrl_1', 'ctrl_2', 'treat_1', 'treat_2'),
    treatment = c(0, 0, 1, 1)
)
```

### 2. Read Data

```r
library(methylKit)

meth_obj <- methRead(
    location = as.list(samples$file),
    sample.id = as.list(samples$sample_id),
    treatment = samples$treatment,
    assembly = 'hg38',
    pipeline = 'bismarkCoverage'
)
```

### 3. Quality Control

```r
# Check each sample
for (i in seq_along(meth_obj)) {
    getMethylationStats(meth_obj[[i]], plot = TRUE)
    getCoverageStats(meth_obj[[i]], plot = TRUE)
}
```

### 4. Filter and Normalize

```r
meth_filt <- filterByCoverage(meth_obj, lo.count = 10, hi.perc = 99.9)
meth_norm <- normalizeCoverage(meth_filt)
meth_united <- unite(meth_norm, destrand = TRUE)
```

### 5. Sample Visualization

```r
getCorrelation(meth_united, plot = TRUE)
PCASamples(meth_united)
```

### 6. Differential Methylation

```r
diff_meth <- calculateDiffMeth(meth_united, overdispersion = 'MN', mc.cores = 4)
dmcs <- getMethylDiff(diff_meth, difference = 25, qvalue = 0.01)
```

## Common Issues

### Low Number of Common CpGs

- Relax coverage filter (lo.count = 5)
- Use min.per.group in unite()
- Ensure samples are from same protocol

### Memory Issues

- Process chromosomes separately
- Use save.db = TRUE for database-backed objects

## Output Interpretation

| Column | Description |
|--------|-------------|
| chr, start, end | Genomic position |
| meth.diff | Methylation difference (%) |
| pvalue | Raw p-value |
| qvalue | FDR-adjusted p-value |

Positive meth.diff = hypermethylated in treatment
Negative meth.diff = hypomethylated in treatment

## Resources

- [methylKit Bioconductor](https://bioconductor.org/packages/methylKit/)
- [methylKit Vignette](https://bioconductor.org/packages/release/bioc/vignettes/methylKit/inst/doc/methylKit.html)

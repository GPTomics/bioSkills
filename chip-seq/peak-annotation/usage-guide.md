# Peak Annotation with ChIPseeker - Usage Guide

## Overview

ChIPseeker annotates ChIP-seq peaks to genomic features and genes. It provides statistical summaries, visualizations, and functional enrichment analysis of peak-associated genes.

## When to Use This Skill

- You have peak files from MACS2 or other peak callers
- You want to know which genes are associated with peaks
- You need to visualize peak distribution across genomic features
- You want to compare peak annotations between samples

## Installation

```r
if (!require('BiocManager', quietly = TRUE))
    install.packages('BiocManager')

BiocManager::install(c('ChIPseeker', 'TxDb.Hsapiens.UCSC.hg38.knownGene', 'org.Hs.eg.db'))

# For functional enrichment
BiocManager::install('clusterProfiler')
```

## Basic Workflow

### 1. Load Packages and Data

```r
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peaks <- readPeakFile('sample_peaks.narrowPeak')
```

### 2. Annotate Peaks

```r
peak_anno <- annotatePeak(peaks, TxDb = txdb, annoDb = 'org.Hs.eg.db')
```

### 3. Visualize

```r
plotAnnoPie(peak_anno)
plotDistToTSS(peak_anno)
```

### 4. Export Results

```r
anno_df <- as.data.frame(peak_anno)
write.csv(anno_df, 'annotated_peaks.csv', row.names = FALSE)
```

## Common Issues

### Chromosome Name Mismatch

UCSC uses "chr1", Ensembl uses "1". Ensure consistency:
```r
# Add chr prefix if needed
seqlevelsStyle(peaks) <- 'UCSC'
```

### No Gene Symbols

Ensure annoDb parameter is set correctly:
```r
annotatePeak(peaks, TxDb = txdb, annoDb = 'org.Hs.eg.db')
```

## Output Columns

| Column | Description |
|--------|-------------|
| seqnames | Chromosome |
| start, end | Peak coordinates |
| annotation | Genomic feature |
| distanceToTSS | Distance to nearest TSS |
| SYMBOL | Gene symbol |
| GENENAME | Gene description |
| ENTREZID | Entrez gene ID |

## Resources

- [ChIPseeker Bioconductor](https://bioconductor.org/packages/ChIPseeker/)
- [ChIPseeker Vignette](https://bioconductor.org/packages/release/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html)

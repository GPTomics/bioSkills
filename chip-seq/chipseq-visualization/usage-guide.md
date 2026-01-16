# ChIP-seq Visualization - Usage Guide

## Overview

This skill covers visualization of ChIP-seq data using deepTools (CLI), ChIPseeker (R), and Gviz (R). Create heatmaps, profile plots, and genome browser-style visualizations.

## Tool Comparison

| Tool | Type | Use Case |
|------|------|----------|
| deepTools | CLI | Fast, large-scale analysis |
| ChIPseeker | R | Peak-centric, TSS profiles |
| Gviz | R | Publication-quality browser |
| EnrichedHeatmap | R | Customizable heatmaps |

## Installation

### deepTools (CLI)

```bash
conda install -c bioconda deeptools
```

### R Packages

```r
BiocManager::install(c('ChIPseeker', 'Gviz', 'EnrichedHeatmap'))
```

## Basic Workflows

### Heatmap with deepTools

```bash
# 1. Create bigWig from BAM
bamCoverage -b sample.bam -o sample.bw --normalizeUsing CPM

# 2. Compute matrix
computeMatrix reference-point -R genes.bed -S sample.bw -o matrix.gz

# 3. Plot heatmap
plotHeatmap -m matrix.gz -o heatmap.png
```

### Profile with ChIPseeker

```r
library(ChIPseeker)
peaks <- readPeakFile('peaks.narrowPeak')
tagMatrix <- getTagMatrix(peaks, windows = promoter)
plotAvgProf(tagMatrix, xlim = c(-3000, 3000))
```

### Browser View with Gviz

```r
library(Gviz)
dtrack <- DataTrack(range = 'sample.bw', genome = 'hg38')
plotTracks(dtrack, from = 1000000, to = 1100000, chromosome = 'chr1')
```

## Common bigWig Normalization

| Method | Description |
|--------|-------------|
| CPM | Counts per million |
| RPKM | Reads per kilobase million |
| BPM | Bins per million |
| RPGC | Reads per genomic content |

## Resources

- [deepTools Documentation](https://deeptools.readthedocs.io/)
- [ChIPseeker Vignette](https://bioconductor.org/packages/ChIPseeker/)
- [Gviz User Guide](https://bioconductor.org/packages/Gviz/)

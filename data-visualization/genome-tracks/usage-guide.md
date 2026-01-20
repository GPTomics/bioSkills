# Genome Track Visualization Usage Guide

## Overview

Genome track plots show multiple data layers (coverage, peaks, genes) aligned to genomic coordinates, similar to genome browsers.

## Tool Selection

| Tool | Language | Best For |
|------|----------|----------|
| pyGenomeTracks | Python/CLI | Publication figures, batch processing |
| Gviz | R | R workflows, Bioconductor integration |
| IGV.js | JavaScript | Web embedding, interactive |

## Quick Start Prompts

- "Create a genome track plot with coverage and peaks"
- "Show gene models with ChIP-seq signal above"
- "Compare two samples side by side at a locus"
- "Plot multiple genomic regions from a BED file"

## Requirements

```bash
# pyGenomeTracks
pip install pygenometracks

# R/Gviz
BiocManager::install('Gviz')
```

## Input Files

- **BigWig** (.bw) - Continuous signal (coverage)
- **BED** (.bed) - Intervals (peaks, regions)
- **GTF/GFF** - Gene annotations
- **narrowPeak** - MACS3 peak files

## Related Skills

- **genome-intervals/bigwig-files** - Create BigWig files
- **chip-seq/chipseq-visualization** - ChIP-specific plots

# Genome Tracks - Usage Guide

## Overview
Genome track plots show multiple data layers (coverage, peaks, genes) aligned to genomic coordinates, similar to genome browsers.

## Prerequisites
```bash
# pyGenomeTracks
pip install pygenometracks

# R/Gviz
BiocManager::install('Gviz')
```

## Quick Start
Tell your AI agent what you want to do:
- "Create a genome track plot with coverage and peaks"
- "Show gene models with ChIP-seq signal above"
- "Compare two samples side by side at a locus"

## Example Prompts
### Basic Tracks
> "Create a track plot showing my BigWig coverage at the TP53 locus"

> "Plot gene models from my GTF file"

### Multi-Track Figures
> "Show ChIP-seq signal with peaks and gene annotations below"

> "Compare treatment vs control coverage at this region"

### Batch Processing
> "Plot multiple genomic regions from my BED file"

> "Create track plots for all genes in my list"

## What the Agent Will Do
1. Set up track configuration (BigWig, BED, GTF files)
2. Define the genomic region to visualize
3. Configure track heights and colors
4. Generate the track plot
5. Export at publication resolution

## Tool Selection

| Tool | Language | Best For |
|------|----------|----------|
| pyGenomeTracks | Python/CLI | Publication figures, batch processing |
| Gviz | R | R workflows, Bioconductor integration |
| IGV.js | JavaScript | Web embedding, interactive |

## Input Files
- **BigWig** (.bw) - Continuous signal (coverage)
- **BED** (.bed) - Intervals (peaks, regions)
- **GTF/GFF** - Gene annotations
- **narrowPeak** - MACS3 peak files

## Tips
- Use BigWig files for efficient coverage visualization
- Match y-axis scales when comparing samples
- Include gene annotations for context
- Use pyGenomeTracks config files for reproducibility

## Related Skills
- **genome-intervals/bigwig-files** - Create BigWig files
- **chip-seq/chipseq-visualization** - ChIP-specific plots

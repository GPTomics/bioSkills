# ChIP-seq Visualization - Usage Guide

## Overview

Create heatmaps, profile plots, and genome browser-style visualizations for ChIP-seq data using deepTools (CLI), ChIPseeker (R), and Gviz (R).

## Prerequisites

```bash
# deepTools (CLI)
conda install -c bioconda deeptools

# R packages
# BiocManager::install(c('ChIPseeker', 'Gviz', 'EnrichedHeatmap'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Create a heatmap of ChIP-seq signal around TSS regions"
- "Generate profile plots comparing treatment vs control"
- "Make a genome browser view of my peaks at the MYC locus"

## Example Prompts

### Heatmaps and Profiles
> "Create a heatmap showing H3K27ac signal at all promoters"

> "Compare ChIP-seq signal profiles between wild-type and knockout samples"

> "Generate a reference-point matrix centered on peak summits"

### Genome Browser Views
> "Create a publication-quality browser track for chromosome 1:1-2Mb"

> "Visualize my bigWig files alongside gene annotations"

### Normalization
> "Convert my BAM files to bigWig with CPM normalization"

> "Create RPGC-normalized coverage tracks for my samples"

## What the Agent Will Do

1. Convert BAM files to bigWig format with appropriate normalization (CPM, RPKM, RPGC)
2. Compute signal matrices around reference points (TSS, peak summits) or scaled regions
3. Generate heatmaps and profile plots using deepTools or EnrichedHeatmap
4. Create genome browser visualizations with Gviz for specific loci
5. Customize plot aesthetics and export publication-ready figures

## Tips

- Use CPM normalization for comparing samples with similar library sizes
- RPGC (reads per genomic content) is better for comparing samples with different sequencing depths
- For TSS profiles, use a window of +/- 3kb around the TSS
- deepTools is faster for large-scale analysis; R packages offer more customization
- Always include a control/input track when visualizing ChIP-seq data

# Circos Plots - Usage Guide

## Overview
Create circular genome visualizations with multiple data tracks for CNVs, gene fusions, Hi-C contacts, and multi-omics summaries.

## Prerequisites
```bash
# Circos (Perl)
conda install -c bioconda circos

# pyCircos (Python)
pip install pyCircos

# circlize (R)
install.packages('circlize')
```

## Quick Start
Tell your AI agent what you want to do:
- "Create a circos plot of my CNV data"
- "Visualize gene fusions as links between chromosomes"
- "Make a circular genome plot with expression data"

## Example Prompts
### CNV Visualization
> "Create a circos plot showing copy number gains in red and losses in blue"

> "Add a CNV heatmap track to my genome circle"

### Gene Fusions
> "Draw arcs connecting fusion breakpoints between chromosomes"

> "Visualize structural variants as links in a circos plot"

### Hi-C and Chromatin
> "Show Hi-C contact frequencies as links between regions"

> "Create a circos plot of chromatin interactions"

### Multi-omics Summary
> "Layer variant, expression, and methylation data around the genome"

> "Create a summary circos with multiple data tracks"

## What the Agent Will Do
1. Initialize genome coordinates (hg38/mm10)
2. Add chromosome ideograms as outer track
3. Add data tracks (bars, heatmaps, scatter)
4. Draw links for interactions/fusions
5. Export as SVG/PNG

## Common Use Cases
- **CNV Visualization**: Show copy number gains (red) and losses (blue) as bars or heatmap tracks
- **Gene Fusions**: Display fusion breakpoints as arcs connecting two genomic positions
- **Hi-C Contacts**: Visualize chromatin interactions as links between genomic regions
- **Multi-omics Summary**: Layer multiple data types around the genome

## Tips
- Track order: Outermost = chromosome ideograms, inner tracks = data
- Use `interspace` parameter to separate chromosomes
- Keep color schemes consistent across tracks
- Export as SVG for publication, PNG for quick viewing
- Use circlize (R) for fastest prototyping

## Related Skills
- **data-visualization/genome-tracks** - Linear track plots
- **hi-c-analysis/contact-visualization** - Hi-C specific plots

# data-visualization

## Overview

Publication-quality data visualization for bioinformatics using ggplot2 and matplotlib with best practices for scientific figures.

**Tool type:** mixed | **Primary tools:** ggplot2, matplotlib, plotly, ComplexHeatmap

## Skills

| Skill | Description |
|-------|-------------|
| ggplot2-fundamentals | Create publication-ready figures with ggplot2 |
| bindplots-customization | Combine and customize multi-panel figures |
| heatmaps-clustering | Expression heatmaps with ComplexHeatmap and pheatmap |
| interactive-visualization | Interactive plots with plotly and bokeh |
| genome-tracks | Genome browser tracks with pyGenomeTracks and Gviz |
| specialized-omics-plots | Volcano, MA, PCA, and enrichment dotplots |
| color-palettes | Colorblind-friendly palettes and journal color schemes |

## Example Prompts

- "Create a publication-quality volcano plot"
- "Make a multi-panel figure with shared legends"
- "Set up a consistent color scheme for my figures"
- "Export figures at 300 DPI for publication"
- "Create an interactive heatmap for my expression data"
- "Plot genome tracks for my ChIP-seq regions"
- "Apply a colorblind-friendly palette"

## Requirements

```r
# R packages
install.packages(c('ggplot2', 'patchwork', 'scales', 'ggrepel', 'pheatmap'))
BiocManager::install(c('ComplexHeatmap', 'Gviz'))
```

```bash
# Python packages
pip install matplotlib seaborn plotly bokeh pyGenomeTracks
```

## Related Skills

- **differential-expression/de-visualization** - Expression-specific plots
- **pathway-analysis/enrichment-visualization** - Enrichment plots
- **reporting** - Figures in reports

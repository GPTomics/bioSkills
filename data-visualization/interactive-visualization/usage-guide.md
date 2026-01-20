# Interactive Visualization Usage Guide

## Overview

Interactive plots enable exploration of large datasets through zooming, panning, and hover tooltips. Export as standalone HTML files for sharing.

## Tool Selection

| Tool | Language | Best For |
|------|----------|----------|
| plotly | Python/R | General interactive plots, ggplot2 conversion |
| bokeh | Python | Web apps, linked brushing, widgets |

## Quick Start Prompts

- "Create an interactive volcano plot with gene hover labels"
- "Make a zoomable PCA plot colored by condition"
- "Build an interactive heatmap of my expression data"
- "Convert my ggplot to interactive plotly"

## Requirements

```bash
# Python
pip install plotly bokeh

# R
install.packages('plotly')
```

## Key Features

- **Hover tooltips** - Show metadata on mouseover
- **Zoom/pan** - Explore dense regions
- **Selection** - Click to select points
- **Export** - Standalone HTML files

## Related Skills

- **data-visualization/ggplot2-fundamentals** - Static versions
- **reporting/quarto-reports** - Embed in documents

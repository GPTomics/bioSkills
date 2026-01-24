# Spatial Visualization - Usage Guide

## Overview

This skill covers creating visualizations for spatial transcriptomics data using Squidpy and Scanpy, including tissue plots with gene expression and cluster overlays.

## Prerequisites

```bash
pip install squidpy scanpy matplotlib
# Optional for interactive:
pip install napari
```

## Quick Start

Tell your AI agent what you want to do:
- "Plot gene expression on the tissue"
- "Show clusters on the spatial plot"

## Example Prompts

### Gene Expression
> "Plot CD3D expression on the tissue"

> "Show expression of these marker genes spatially"

### Clusters
> "Visualize clusters on the tissue image"

> "Color spots by cluster"

### Publication Figures
> "Create a multi-panel figure with clusters and markers"

> "Save a high-resolution spatial plot"

## What the Agent Will Do

1. Load spatial data with coordinates
2. Create matplotlib figure
3. Plot spots colored by requested variable
4. Optionally overlay on tissue image
5. Save or display figure

## Tips

- **spot_size** - Adjust based on spot density
- **vmin/vmax** - Use percentiles (e.g., 'p99') for outlier-robust scaling
- **alpha_img** - Control tissue image transparency
- **frameon=False** - Cleaner look for publication figures

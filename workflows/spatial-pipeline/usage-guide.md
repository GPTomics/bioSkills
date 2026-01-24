# Spatial Transcriptomics Pipeline - Usage Guide

## Overview

This workflow analyzes spatial transcriptomics data (Visium, Xenium) from raw data to spatial domains and visualizations using Squidpy and Scanpy.

## Prerequisites

```bash
pip install squidpy scanpy matplotlib
```

## Quick Start

Tell your AI agent what you want to do:
- "Analyze my Visium spatial transcriptomics data"
- "Find spatially variable genes in my tissue"
- "Identify spatial domains in my sample"

## Example Prompts

### Loading and QC
> "Load my Space Ranger output"

> "Show QC metrics on the tissue image"

### Analysis
> "Find spatially variable genes"

> "Run neighborhood enrichment analysis"

> "Detect spatial domains"

### Visualization
> "Plot gene expression on the tissue"

> "Show clusters overlaid on the image"

## Input Requirements

| Input | Format | Description |
|-------|--------|-------------|
| Space Ranger output | Directory | Visium processed data |
| Xenium output | Directory | Xenium processed data |

## What the Workflow Does

1. **Load Data** - Read spatial data with images
2. **QC** - Filter low-quality spots
3. **Clustering** - Standard scRNA-seq pipeline
4. **Spatial Analysis** - Neighbors, statistics
5. **Domains** - Spatial domain detection
6. **Visualization** - Plots on tissue

## Tips

- **Spot size**: Adjust for visualization clarity
- **Resolution**: Lower for fewer, larger domains
- **SVGs**: Check top Moran's I genes
- **Deconvolution**: Add for cell type estimates

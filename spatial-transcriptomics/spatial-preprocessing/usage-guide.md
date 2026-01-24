# Spatial Preprocessing - Usage Guide

## Overview

This skill covers quality control, filtering, normalization, and feature selection for spatial transcriptomics data using Squidpy and Scanpy.

## Prerequisites

```bash
pip install squidpy scanpy matplotlib
```

## Quick Start

Tell your AI agent what you want to do:
- "Run QC on my spatial data"
- "Filter and normalize my Visium data"

## Example Prompts

### QC
> "Calculate QC metrics for my spatial data"

> "Show QC metrics on the tissue"

### Filtering
> "Filter spots with less than 500 counts"

> "Remove spots with high mitochondrial content"

### Normalization
> "Normalize my spatial data"

> "Find highly variable genes"

## What the Agent Will Do

1. Calculate QC metrics (counts, genes, MT%)
2. Visualize QC metrics on tissue
3. Filter low-quality spots
4. Normalize expression data
5. Identify highly variable genes
6. Optionally find spatially variable genes

## Tips

- **Spatial QC** - Always visualize QC metrics on the tissue to identify spatial artifacts
- **Mitochondrial threshold** - Often higher for spatial data (~20-25%)
- **SVG vs HVG** - Spatially variable genes may differ from highly variable genes
- **Keep raw counts** - Store in `adata.layers['counts']` before normalization

# Spatial Transcriptomics Pipeline - Usage Guide

## Overview

This workflow analyzes spatial transcriptomics data (Visium, Xenium) from raw data to spatial domains and visualizations using Squidpy and Scanpy. It branches first on the platform-class fork: sequencing/capture data (Visium) are spot mixtures that are QC'd, clustered into regions, and deconvolved for composition; imaging data (Xenium) are single cells that need low-count QC floors and segmentation. Composition estimation (spatial-deconvolution) and cell-cell communication (spatial-communication) are handed off to those skills rather than inlined.

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

1. **Load Data** - Read spatial data with images (squidpy for Visium; spatialdata_io for Xenium, which has no squidpy reader)
2. **QC** - Filter spots with platform-appropriate floors (Visium floors delete real imaging cells, so branch on the fork)
3. **Clustering** - Cluster spots into regions/niches, not cell types (a Visium spot is a 1-10-cell mixture)
4. **Spatial Analysis** - Neighbor graph on the hex lattice (coord_type grid), Moran's I gated on FDR, neighborhood enrichment read against its weak null
5. **Domains** - Spatial domain detection (a domain is a region, distinct from a cell type and a niche)
6. **Visualization** - Plots on tissue with the correct spot size and coordinate frame

## Tips

- **Platform fork**: decide imaging vs sequencing first; it sets QC floors and whether you deconvolve or segment
- **QC floors**: Visium UMI/gene floors are tissue-dependent starting points; never apply a 500-count floor to imaging data
- **SVGs**: gate top Moran's I on FDR, and separate composition-driven SVGs (cell-type markers) from within-type regulation
- **nhood enrichment**: a positive z is not a specific interaction; the global null is confounded by abundance and compartments
- **Deconvolution and communication**: run them through spatial-deconvolution and spatial-communication, which carry the reference and spillover caveats

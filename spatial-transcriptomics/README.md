# spatial-transcriptomics

## Overview

Analyze spatial transcriptomics data from Visium, Xenium, MERFISH, and other platforms using Squidpy and SpatialData.

**Tool type:** python | **Primary tools:** Squidpy, SpatialData, Scanpy

## Skills

| Skill | Description |
|-------|-------------|
| spatial-data-io | Load spatial data from Visium, Xenium, Slide-seq, MERFISH |
| spatial-preprocessing | QC, normalization, and feature selection for spatial data |
| spatial-neighbors | Build spatial neighbor graphs and compute connectivity |
| spatial-statistics | Moran's I, spatial autocorrelation, co-occurrence, enrichment |
| spatial-domains | Identify spatial domains and tissue regions |
| image-analysis | Process and analyze tissue images with Squidpy |
| spatial-visualization | Static and interactive visualization of spatial data |
| spatial-communication | Ligand-receptor analysis and cell-cell interactions |
| spatial-deconvolution | Estimate cell type composition per spot |
| spatial-multiomics | Analyze high-resolution platforms (Slide-seq, Stereo-seq, Visium HD) |

## Example Prompts

- "Load my Visium data"
- "Read this Xenium output folder"
- "Run QC on my spatial data"
- "Normalize my spatial transcriptomics data"
- "Build a spatial neighbor graph with 6 neighbors"
- "Calculate Moran's I for this gene"
- "Find spatially variable genes"
- "Run co-occurrence analysis"
- "Identify spatial domains in my tissue"
- "Segment cells from the H&E image"
- "Plot gene expression on the tissue"
- "Show clusters overlaid on the image"
- "Run ligand-receptor analysis"
- "Deconvolve my Visium data with cell2location"
- "Analyze my Slide-seq data"
- "Process Stereo-seq at bin level"
- "Work with Visium HD subcellular resolution"

## Requirements

```bash
pip install squidpy spatialdata spatialdata-io scanpy anndata
```

## Related Skills

- **single-cell** - Non-spatial scRNA-seq analysis
- **image-analysis** - General image processing
- **differential-expression** - DE between spatial regions

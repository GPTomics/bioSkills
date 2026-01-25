# spatial-transcriptomics

## Overview

Analyze spatial transcriptomics data from Visium, Xenium, MERFISH, and other platforms using Squidpy and SpatialData.

**Tool type:** python | **Primary tools:** Squidpy, SpatialData, Scanpy, scimap

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
| spatial-proteomics | Analyze CODEX, IMC, MIBI spatial proteomics data |

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
- "Analyze my CODEX spatial proteomics data"
- "Find spatial interactions between cell types in IMC data"

## Requirements

```bash
pip install squidpy spatialdata spatialdata-io scanpy anndata scimap
```

## Related Skills

- **single-cell** - Non-spatial scRNA-seq analysis
- **differential-expression** - DE between spatial regions
- **data-visualization** - Visualization of spatial patterns

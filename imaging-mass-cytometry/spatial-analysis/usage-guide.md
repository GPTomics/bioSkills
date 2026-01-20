# Spatial Analysis - Usage Guide

## Overview
Spatial analysis reveals tissue organization by analyzing cell type distributions and interactions in their spatial context.

## Prerequisites
```bash
pip install squidpy scanpy anndata scipy
```

## Quick Start
Tell your AI agent what you want to do:
- "Analyze spatial relationships between cell types"
- "Test if T cells cluster around tumor cells"
- "Calculate neighborhood enrichment scores"

## Example Prompts

### Neighborhood Enrichment
> "Calculate neighborhood enrichment between all cell types in my tissue"

> "Test if CD8 T cells are enriched near tumor cells"

> "Run permutation-based neighborhood analysis on my IMC data"

### Co-occurrence Analysis
> "Analyze how cell type co-occurrence changes with distance"

> "Calculate co-occurrence curves for immune and tumor cells"

### Spatial Clustering
> "Test if macrophages are clustered or randomly distributed"

> "Calculate Ripley's L function for my cell populations"

### Spatial Graph Construction
> "Build a spatial neighborhood graph using Delaunay triangulation"

> "Create a radius-based spatial graph with 50 micron cutoff"

## What the Agent Will Do
1. Load cell coordinates and phenotype annotations
2. Build spatial graph (Delaunay or radius-based)
3. Calculate neighborhood enrichment z-scores via permutation
4. Optionally compute co-occurrence across distance ranges
5. Generate spatial statistics plots (enrichment heatmaps, co-occurrence curves)
6. Report significant cell-cell interactions

## Tips
- Delaunay triangulation works well for densely packed tissues
- Radius-based graphs (2-3 cell diameters) are better for sparse tissues
- Z-score > 2: significant attraction; Z-score < -2: significant avoidance
- Always normalize for cell type abundance when interpreting enrichment
- Test multiple radii to assess robustness of spatial patterns
- squidpy provides comprehensive spatial analysis functions for AnnData objects

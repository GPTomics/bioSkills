# Spatial Proteomics Usage Guide

## Overview

This guide covers analyzing spatial proteomics data from CODEX, IMC, and MIBI platforms.

## Prerequisites

```bash
pip install scimap squidpy anndata scanpy
```

## Quick Start

Tell your AI agent what you want to do:
- "Load and preprocess my CODEX data"
- "Phenotype cells based on marker expression"
- "Analyze spatial interactions between cell types"
- "Find cellular neighborhoods in my IMC data"

## Example Prompts

### Data Processing

> "Load my spatial proteomics h5ad file and normalize the marker intensities"

> "Apply batch correction across my multiple FOVs"

### Cell Phenotyping

> "Phenotype cells as T cells, B cells, macrophages, and tumor cells based on canonical markers"

> "Cluster my cells and annotate clusters based on marker expression"

### Spatial Analysis

> "Calculate spatial interactions between tumor and immune cells"

> "Identify cellular neighborhoods and their composition"

> "Find which cell types are enriched near tumor boundaries"

### Visualization

> "Create a spatial plot colored by cell phenotype"

> "Generate a heatmap of cell-cell spatial interactions"

## What the Agent Will Do

1. Load spatial proteomics data into AnnData format
2. Preprocess and normalize marker intensities
3. Perform cell phenotyping (gating or clustering)
4. Build spatial neighbor graphs
5. Calculate spatial statistics and interactions
6. Generate visualizations

## Tips

- Normalize markers before phenotyping
- Use rescale to make markers comparable
- Set reasonable gating thresholds (often 0.5 after rescale)
- Consider k-nearest neighbors (k=10-30) for spatial analysis
- Compare spatial patterns across multiple FOVs
- Integrate with transcriptomics if available

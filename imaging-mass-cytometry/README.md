# imaging-mass-cytometry

## Overview

Spatial proteomics analysis from imaging mass cytometry (IMC) and multiplexed ion beam imaging (MIBI) data.

**Tool type:** python | **Primary tools:** steinbock, squidpy, napari

## Skills

| Skill | Description |
|-------|-------------|
| data-preprocessing | Load and preprocess IMC/MIBI data |
| cell-segmentation | Segment cells from multiplexed images |
| phenotyping | Assign cell types from marker expression |
| spatial-analysis | Analyze spatial relationships and neighborhoods |
| interactive-annotation | Manual cell type annotation with napari |
| quality-metrics | QC metrics for SNR, drift, and artifacts |

## Example Prompts

- "Load my IMC data and segment cells"
- "Phenotype cells based on marker expression"
- "Analyze spatial neighborhoods in my tissue"
- "Find spatial interactions between cell types"

## Requirements

```bash
# Python
pip install steinbock squidpy napari scikit-image readimc tifffile anndata scanpy

# Docker (for steinbock)
docker pull ghcr.io/bodenmillergroup/steinbock:latest
```

## Related Skills

- **spatial-transcriptomics** - Similar spatial analysis concepts
- **single-cell** - Cell type annotation methods
- **flow-cytometry** - Similar marker-based phenotyping

# imaging-mass-cytometry

## Overview

Spatial proteomics analysis from imaging mass cytometry (IMC) and multiplexed ion beam imaging (MIBI) data, treating pixels as the ion counts they are, segmentation as the dominant error source, and the patient as the experimental unit.

**Tool type:** mixed | **Primary tools:** steinbock, CATALYST, deepcell, squidpy, napari, diffcyt

## Skills

| Skill | Description |
|-------|-------------|
| data-preprocessing | Ingest MCD/TXT, remove hot pixels, NNLS spillover compensation, arcsinh-cofactor-1 transformation |
| cell-segmentation | Whole-cell/nuclear segmentation (Mesmer, Cellpose, ilastik), error propagation, lateral spillover |
| phenotyping | Assign cell types via clustering or marker/image classifiers; the double-positive artifact |
| spatial-analysis | Neighborhood, niche, and interaction analysis with explicit nulls and graph choices |
| differential-analysis | Compare composition and spatial features across conditions with the patient as the unit |
| interactive-annotation | napari/Mantis annotation and segmentation QC; the pixels-to-cell-table bridge |
| quality-metrics | Multi-level QC: Poisson-count SNR, spillover-matrix QC, drift, and batch effects |

## Example Prompts

- "Load my MCD files, compensate channel spillover with NNLS, and segment cells"
- "Phenotype my cells and tell me if my CD3+CD20+ population is a segmentation artifact"
- "Test whether a niche is real or just a density gradient"
- "Compare cell-type proportions between conditions at the patient level"
- "Decide which channels and ROIs to drop before analysis"

## Requirements

```bash
# Python
pip install steinbock readimc deepcell cellpose squidpy scanpy anndata napari napari-imc statsmodels sccoda numpy

# R / Bioconductor
# BiocManager::install(c('CATALYST', 'cytomapper', 'imcRtools', 'diffcyt', 'spillR'))

# Docker (steinbock)
docker pull ghcr.io/bodenmillergroup/steinbock
```

## Related Skills

- **spatial-transcriptomics** - shared squidpy neighborhood and spatial-statistics methods
- **single-cell** - cell-type annotation and clustering concepts
- **flow-cytometry** - suspension spillover, FlowSOM phenotyping, and diffcyt differential testing
- **experimental-design** - experimental-unit and pseudoreplication foundations

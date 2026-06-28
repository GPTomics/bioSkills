# Spatial Proteomics Usage Guide

## Overview

This guide covers analysis of multiplexed antibody-imaging data -- CODEX/PhenoCycler, MIBI-TOF, IMC, CyCIF, and Opal/Vectra mIF. The central idea is that the signal is continuous protein intensity with antibody, batch, and staining confounds, not transcript counts, so the transform (arcsinh, not log1p-of-counts), channel-spillover compensation, antibody-batch normalization, and intensity-based phenotyping all differ from RNA-based spatial analysis. The antibody panel is targeted (absence is uninformative) and whole-cell segmentation is the dominant downstream error source. For the IMC/MIBI end-to-end processing pipeline, this skill cross-references the imaging-mass-cytometry category rather than duplicating it.

## Prerequisites

```bash
pip install scimap squidpy scanpy anndata
```

## Quick Start

Tell your AI agent what you want to do:
- "Transform my CODEX marker intensities with arcsinh and correct antibody batch"
- "Phenotype my cells by gating on canonical lineage markers"
- "Cluster my multiplexed-imaging cells and annotate against the panel"
- "Test which cell types are spatial neighbors more than chance"
- "Find recurrent cellular neighborhoods in my MIBI data"

## Example Prompts

### Transform and Normalize
> "My data is CODEX protein intensity -- pick and apply a variance-stabilizing transform and tune the arcsinh cofactor so near-zero noise is not over-expanded."

> "Correct antibody-batch and staining-day differences across my images before I compare samples."

### Cell Phenotyping
> "Phenotype cells as T cells, B cells, macrophages, and tumor using gating on my lineage markers, and flag any double-positive cells that look like segmentation artifacts."

> "Cluster my cells on transformed intensities and annotate the clusters, keeping in mind the panel is targeted."

### Spatial Analysis
> "Build a spatial graph on cell centroids and run a permutation neighborhood-enrichment test between tumor and immune cells."

> "Summarize each cell's local window into cellular neighborhoods and report sensitivity to the window size."

### Platform and Boundary
> "My data is IMC -- decide whether to handle spillover compensation and segmentation here or in the imaging-mass-cytometry pipeline."

## What the Agent Will Do

1. Identify the platform (CODEX/MIBI/IMC/CyCIF/Opal) and its dominant confounds.
2. Apply an intensity transform (arcsinh with a tuned cofactor, or z-score/percentile) -- never log1p-of-counts.
3. Compensate channel spillover and correct antibody-batch/staining effects before comparison.
4. Phenotype cells by gating (workflow DataFrame) or clustering on transformed intensities.
5. Build a spatial graph and run permutation-based neighborhood-enrichment and niche analysis.
6. Audit segmentation-driven artifacts (phantom double-positives) and panel-bounded absences.

## Tips

- Treat intensity as continuous and confounded: arcsinh (tune the cofactor; ~5 is a CyTOF convention, not auto-optimal for imaging), z-score, or percentile -- never log1p applied as if values were UMIs.
- Channel spillover and segmentation lateral spillover both fabricate double-positive cells; compensate and audit segmentation before trusting any phenotype.
- The antibody panel is targeted: a marker or cell type that is "absent" was most likely not stained, so absence is uninformative -- there is no de-novo discovery.
- Normalization x clustering choices dominate cell-type calls (Hickey 2021: one CODEX dataset gave 20 annotations); fix and report the transform, cofactor, normalization, and clustering.
- Whole-cell segmentation on a membrane stain (Mesmer) is the highest-leverage decision; every neighborhood result inherits segmentation error.
- For the IMC/MIBI end-to-end pipeline (NNLS spillover, segmentation execution, FlowSOM phenotyping on metal channels), use the imaging-mass-cytometry skills.

## Related Skills

- image-analysis - whole-cell segmentation upstream of every per-cell intensity vector
- imaging-mass-cytometry/cell-segmentation - segmentation execution and error propagation for IMC/MIBI
- imaging-mass-cytometry/phenotyping - FlowSOM/Phenograph phenotyping on metal channels
- spatial-transcriptomics/spatial-multiomics - integrating spatial proteomics with matched transcriptomics
- spatial-transcriptomics/spatial-statistics - permutation nulls and neighborhood enrichment shared with squidpy

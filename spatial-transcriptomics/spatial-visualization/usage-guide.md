# Spatial Visualization - Usage Guide

## Overview

This skill covers plotting spatial transcriptomics expression, clusters, and continuous scores on a tissue section with Squidpy and Scanpy. The central decisions are which plotter and marker size to use for the platform fork (spot/capture like Visium and Slide-seq versus imaging/FOV like Xenium, MERFISH, and CosMx), how to align points to the histology image with the correct coordinate-frame transform, and how to avoid rendering choices -- interpolation, oversized markers, misleading colormaps, non-metric embedding distances -- that manufacture structure the assay never measured.

## Prerequisites

```bash
pip install squidpy scanpy anndata matplotlib
# Optional perceptually-uniform colormaps and interactive viewers:
pip install cmcrameri napari vitessce
```

## Quick Start

Tell your AI agent what you want to do:
- "Plot CD3D expression on the Visium tissue with the histology image"
- "Show clusters on my Xenium section as a cell point cloud"
- "Color spots by a continuous score without smoothing the field"
- "Make a publication panel with a perceptually-uniform colormap"

## Example Prompts

### Spot/Capture (Visium, Slide-seq)
> "Overlay leiden clusters on the H&E image using the dataset scalefactors and a spot size matching the capture diameter"

> "Plot total counts and these marker genes side by side on the tissue"

### Imaging/FOV (Xenium, MERFISH, CosMx)
> "Plot cell types as a point cloud in micron coordinates without a spot lattice"

> "Overlay transcripts on the DAPI image with the correct micron-to-pixel transform"

### Honest Continuous Fields
> "Show this spatial score on the measured spots only, no interpolation, with a perceptually-uniform colormap and a disclosed p99 clip"

### Interactive
> "Open this section in napari with the tissue image and spots aligned in pixel space"

## What the Agent Will Do

1. Determine the platform fork (spot/capture vs imaging/FOV) and pick the matching plotter and marker-size meaning.
2. Resolve the coordinate-frame transform (scalefactors for spot/capture, the platform affine for imaging) before any histology overlay.
3. Render measured spots/cells directly, refusing to interpolate, KDE, or contour a sparse field.
4. Choose a perceptually uniform colormap and disclose any vmin/vmax clipping.
5. Restrict spatial conclusions to the spatial plot and treat UMAP/tSNE as QC only.

## Tips

- A capture spot is a 1-10-cell mixture, not a cell -- label spot clusters as regions/niches, not cell types.
- In `scanpy.pl.spatial`, `size` scales the spot diameter; oversizing it fakes tissue coverage the assay did not resolve.
- In `squidpy.pl.spatial_scatter` with no image, `size` is the actual dot size -- keep it small for sparse imaging detections.
- Never KDE/kriging/contour a sparse spatial field; the gradients can be artifacts of the kernel, and spatial statistics on the smoothed surface are partly circular.
- Avoid jet/rainbow colormaps; use viridis or Crameri scientific colour maps (cmcrameri) and disclose any color-scale clipping.
- Always check the histology overlay alignment; a small coordinate-frame error puts expression in the wrong structure.
- UMAP/tSNE distances are not metric -- do not read spatial conclusions from embedding gaps.

## Related Skills

- spatial-data-io - load the platform data and the histology image plus scalefactors that plotting depends on
- spatial-domains - produce the region labels rendered on the section
- spatial-statistics - compute Moran's I / neighborhood enrichment whose results are plotted here
- data-visualization/heatmaps-clustering - general perceptually-uniform colormap and figure conventions
- single-cell/clustering - the non-metric UMAP/tSNE distance caveat that applies to embeddings

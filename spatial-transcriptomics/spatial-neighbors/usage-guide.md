# Spatial Neighbor Graphs - Usage Guide

## Overview

This skill builds the spatial neighbor graph -- the sparse weights matrix W over cell or spot coordinates -- that every downstream spatial statistic inherits. Moran's I, neighborhood enrichment, co-occurrence, and graph-based spatial-domain methods are all functions of W, so the graph is not preprocessing; it is an argument to the statistic. The skill covers choosing the graph type by platform (Visium hex grid vs kNN vs Delaunay vs fixed-radius), handling variable cell density, getting coordinate units right, pruning Delaunay long edges, running the graph sensitivity analysis, and knowing when planar neighbors misrepresent a 3D tissue.

## Prerequisites

```bash
pip install squidpy scanpy anndata numpy
```

## Quick Start

Tell your AI agent what you want to do:
- "Build a spatial neighbor graph for my Visium data"
- "Build a Delaunay contact graph for my Xenium cells and prune the long edges"
- "Build a kNN spatial graph and check whether my Moran's I result is stable across k"
- "Check whether my coordinates are in microns or pixels before I set a radius"

## Example Prompts

### Choosing the graph
> "My data is single-cell-resolution MERFISH. Should I use kNN, Delaunay, or a fixed radius for the neighbor graph, and why does it matter?"

> "Build a hexagonal-lattice neighbor graph for my Visium section using the grid geometry."

### Density and units
> "My tissue has dense lymphoid follicles next to sparse stroma. Will a fixed-radius graph bias my neighborhood enrichment, and what should I use instead?"

> "Confirm whether my spatial coordinates are in microns before I set a 30-micron interaction radius."

### Robustness
> "Rebuild my spatial graph under k=6, 15, 30 and Delaunay, and tell me which of my spatially variable genes are graph-robust."

> "My Delaunay graph has edges crossing a necrotic hole. Prune them by maximum edge length."

## What the Agent Will Do

1. Identify the platform class (capture/grid vs imaging point cloud) to set the default geometry.
2. Confirm the coordinate unit (microns vs pixels) before any distance parameter is used.
3. Build the graph with the matching `coord_type` and parameters, storing connectivities and distances under named keys.
4. Prune Delaunay long edges across tissue gaps where relevant.
5. Inspect degree distribution and connected components for over-density or fragmentation.
6. Optionally rebuild under several adjacency definitions and compare the downstream statistic to flag graph-fragile results.

## Tips

- The graph is the model: Moran's I, neighborhood enrichment, and co-occurrence are all `f(expression, W)`, so the graph choice silently sets the answer -- report results under at least two graphs.
- Squidpy's `nhood_enrichment` builds a Delaunay graph by default; many published z-scores are Delaunay-specific and shift under kNN.
- Fixed-radius graphs bias toward dense regions (more neighbors of everything = fake enrichment); kNN fixes the count but distorts physical distance. There is no free lunch -- pick the distortion the tissue can tolerate.
- Coordinate units decide the graph: a radius in pixels when coords are in microns (or vice versa) silently builds the wrong graph. Visium array row/col is a lattice index, not microns.
- Unpruned Delaunay leaps across necrotic holes and folds; prune by a max edge length (a few cell diameters).
- A Visium spot is a 1-10-cell mixture, not a cell -- spot adjacency is not cell-cell contact.
- A section is one ~5-10 micron plane of a 3D tissue; planar neighbors are not guaranteed 3D neighbors, and truncated cells carry partial profiles.

## Related Skills

- spatial-transcriptomics/spatial-statistics - the neighbor graph is the W fed to Moran's I, neighborhood enrichment, and co-occurrence
- spatial-transcriptomics/spatial-domains - graph-based domain methods inherit this adjacency and its over-smoothing tradeoff
- spatial-transcriptomics/spatial-communication - ligand-receptor proximity tests run on this graph and inherit its density bias
- spatial-transcriptomics/spatial-data-io - load coordinates and confirm their unit before building any graph
- single-cell/clustering - expression-space kNN graphs, the non-spatial counterpart

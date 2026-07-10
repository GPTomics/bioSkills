# Spatial Domain Detection - Usage Guide

## Overview

This skill identifies spatial domains -- contiguous tissue regions of homogeneous expression such as cortical layers, tumor core, or stromal bands -- by combining transcriptional similarity with spatial proximity. It draws the line between a domain (a region containing many cell types), a cell type (one cell's identity), and a niche (the local cell-type composition), then guides choosing a domain method by tissue geometry, tuning the spatial-weight knob to avoid over-smoothing, and choosing the number of domains k as a biological decision.

## Prerequisites

```bash
pip install squidpy scanpy scikit-learn
# optional domain methods:
pip install banksy_py            # neighbor-augmented features
pip install STAGATE_pyG          # graph attention autoencoder (or the TF build)
# BayesSpace / BASS / Banksy are R/Bioconductor:
# BiocManager::install(c('BayesSpace', 'Banksy'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Identify spatial domains in my Visium section"
- "Segment my tissue into regions, not cell types"
- "My domains came out as smooth blobs -- fix the over-smoothing"
- "Pick the number of domains and show me the k+-1 sensitivity"

## Example Prompts

### Domain vs Niche vs Cell Type
> "Is this a domain question or a niche question? I want to know which tissue region each spot belongs to."

> "Separate the cortical layers in my DLPFC section into spatially coherent regions."

### Method Choice by Geometry
> "Which domain method fits continuous laminar cortex versus a high-resolution Xenium tumor section?"

> "My tumor nests are scattered, not contiguous -- which method or parameter handles non-contiguous domains?"

### Tuning the Spatial-Weight Knob
> "My domains are salt-and-pepper -- raise the spatial weight."

> "My regions merged into blobs and the boundaries disappeared -- lower the over-smoothing."

> "Sweep the BANKSY lambda and show me where the boundaries are sharpest."

### Choosing k
> "Set the number of domains for my section and report sensitivity to k plus or minus one."

> "Run BayesSpace with q=7 for cortex and import the labels."

## What the Agent Will Do

1. Establish whether the question is a domain (region), a niche (co-occurring types), or a cell type, and route accordingly.
2. Build a spatial neighbor graph appropriate to the platform (hex grid for Visium, kNN/Delaunay for imaging).
3. Run an expression-only baseline to show the salt-and-pepper failure a domain method corrects.
4. Run a neighbor-augmented or graph-based domain method, exposing the spatial-weight knob.
5. Sweep the spatial weight and inspect boundaries against histology to find the over-smoothing edge.
6. Choose k biologically and report k+-1 sensitivity.
7. Name domains from markers (for ranking/labeling, not inference) and validate against known architecture.

## Tips

- A domain is a REGION with many cell types; a niche is which types co-occur; a cell type is one cell. If the question is "which types co-occur," use neighborhood enrichment in spatial-statistics, not domain segmentation.
- The knob that decides the outcome is the spatial weight (BANKSY lambda, BayesSpace smoothing, SpaGCN histology weight, GNN graph radius), not the clustering algorithm. Too much erases boundaries into blobs; too little reverts to salt-and-pepper.
- k is a biological choice, not a statistic to optimize. Report it and show k+-1.
- Non-contiguous domains (scattered tumor nests) break most methods because the spatial prior assumes contiguity (Yuan 2024). Lower the spatial weight or switch to niche analysis.
- DLPFC is the standard benchmark but it is continuous and laminar, which flatters smoothing methods. Do not over-generalize benchmark rankings to non-laminar tissue.
- No universal winner: GNN methods (STAGATE, GraphST, SEDR) and BayesSpace are robust on continuous Visium; BANKSY/BASS/stLearn lead on high-resolution imaging.
- In multi-sample data, use a multi-sample-aware method (BASS, PRECAST, GraphST integration) or batch-integrate first, or domains will track samples.

## Related Skills

- spatial-neighbors - Build and tune the spatial graph every domain method inherits
- spatial-statistics - Niche/neighborhood enrichment when the question is which cell types co-occur
- spatial-deconvolution - Per-spot cell-type composition, a different question from regions
- single-cell/clustering - Non-spatial clustering and the double-dipping caveat on marker tests

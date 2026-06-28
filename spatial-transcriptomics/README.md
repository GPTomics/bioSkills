# spatial-transcriptomics

## Overview

Analyze spatial transcriptomics and multiplexed-imaging data while respecting the one fork that decides the whole pipeline: the platform class. Imaging / in-situ platforms (Xenium, MERFISH/MERSCOPE, CosMx, seqFISH) read single molecules at subcellular resolution over a TARGETED gene panel, so cells must be SEGMENTED, panels are bounded (absence is uninformative), and deconvolution does not apply. Sequencing / capture platforms (Visium, Visium HD, Slide-seq, Stereo-seq, GeoMx) barcode transcripts to spots or bins that are MIXTURES of cells, so the whole transcriptome is recovered but composition must be DECONVOLVED (or reconstructed up from sub-cellular bins). The first question of any spatial dataset is which side of that fork it sits on; every skill below leads with the branch that applies.

**Tool type:** python | **Primary tools:** Squidpy, SpatialData, cell2location, scimap

## Skills

| Skill | Description |
|-------|-------------|
| spatial-data-io | Load Visium/Visium HD/Xenium/MERFISH/CosMx/Slide-seq/Stereo-seq; pick the reader per platform and separate the per-transcript molecule table from the segmentation-derived cell matrix |
| spatial-preprocessing | Fork-aware QC and normalization; scRNA floors delete real imaging cells, and library size carries spatial biology rather than pure technical depth |
| image-analysis | Segment cells/nuclei (Cellpose/StarDist/Baysor/proseg) and extract image features; judge segmentation spillover before trusting the cell-by-gene matrix |
| high-resolution-binning | Reconstruct single cells from sub-cellular bins (Visium HD 2um, Stereo-seq, Slide-seqV2) with Bin2cell -- bin UP, do not deconvolve DOWN |
| spatial-deconvolution | Estimate per-spot cell-type composition (cell2location/RCTD/SPOTlight); resolution fork (deconvolve vs segment) and reference matching where the reference IS the result |
| spatial-neighbors | Build the spatial neighbor graph (kNN/Delaunay/fixed-radius/Visium hex grid) that every downstream spatial statistic inherits; handle density, units, and the 2D-vs-3D caveat |
| spatial-statistics | Detect spatially variable genes (Moran/Geary/SPARK-X/nnSVG), hot/cold spots (Getis-Ord, LISA), and colocalization with explicit nulls -- separating composition-driven from within-type variation |
| spatial-domains | Identify spatially coherent tissue domains (BANKSY/BayesSpace/STAGATE/GraphST); distinguish a domain from a cell type and a niche, tune the spatial-weight knob, and choose k biologically |
| spatial-communication | Map ligand-receptor cell-cell communication in space; choose a method by whether distance is modeled, pick the L-R database knowingly, and guard the segmentation-spillover circularity |
| spatial-multiomics | Integrate spatial RNA with a second modality (protein/ATAC/histone marks); decide vertical same-pixel (MOFA/WNN) vs diagonal serial-section (PASTE/STalign) integration |
| spatial-proteomics | Analyze multiplexed antibody imaging (CODEX/MIBI/IMC/CyCIF) as continuous protein intensity -- arcsinh transform, spillover/batch correction, gating-vs-clustering phenotyping |
| spatial-visualization | Plot expression, clusters, and scores on tissue with fork-aware plotters and honest rendering (no interpolation, true spot size, perceptually-uniform colormaps) |

## Example Prompts

- "Load my Xenium output and keep both the transcript molecule table and the cell-by-gene matrix"
- "Set QC floors for my Xenium data without deleting real low-count cells"
- "Segment cells from my Xenium data and check whether transcript spillover is faking co-expression"
- "Turn my Visium HD 2um bins into single cells with Bin2cell"
- "Should I deconvolve my Xenium data, or segment it?"
- "Deconvolve my Visium spots into cell-type proportions with cell2location"
- "Build a spatial neighbor graph and check whether my Moran's I result is stable across graph choices"
- "Find spatially variable genes and tell me which are just cell-type markers versus regulated within a cell type"
- "Segment my Visium section into spatial domains and show me the k+-1 sensitivity"
- "Map cell-cell communication and tell me which ligand-receptor calls are spatially supported versus segmentation artifacts"
- "Integrate my spatial CITE-seq RNA and protein into joint factors -- same pixels or adjacent sections?"
- "Analyze my CODEX/MIBI data: transform the protein intensities, phenotype the cells, and find spatial neighbors"
- "Overlay my Xenium cell types as a point cloud and my Visium clusters on the H&E with correct scalefactors"

## Requirements

```bash
pip install squidpy spatialdata spatialdata-io scanpy anndata
pip install cell2location scimap bin2cell
pip install cellpose stardist   # imaging segmentation
```

R-based methods (RCTD/spacexr, SPOTlight, BayesSpace) install via Bioconductor; deep-learning domain methods (STAGATE, GraphST) and Baysor install separately per their docs.

## Related Skills

- **single-cell** - Non-spatial scRNA-seq analysis; the reference for deconvolution and label transfer
- **imaging-mass-cytometry** - The end-to-end IMC/MIBI pipeline (segmentation, spillover, phenotyping); spatial-proteomics covers platform breadth and cross-references it
- **differential-expression** - DE between spatial regions and domains
- **data-visualization** - General plotting patterns reused by spatial-visualization

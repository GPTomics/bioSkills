# Image Analysis - Usage Guide

## Overview

This skill covers cell/nucleus segmentation and image-feature extraction for imaging spatial transcriptomics (Xenium, MERFISH/MERSCOPE, CosMx) and paired H&E/IF tissue images. The central decision is the segmentation strategy: in imaging data there is no native cell -- a cell-by-gene matrix exists only after an algorithm draws boundaries, so segmentation is the dominant upstream error source. The skill helps choose a strategy given the available stain (DAPI nucleus + expansion, membrane-stain whole-cell, transcript-aware Baysor/proseg, or segmentation-free), QC the result for over/under-segmentation and spillover, and decide whether the derived matrix is trustworthy before downstream typing, DE, or ligand-receptor analysis. It also covers Squidpy per-spot image-feature extraction, which is a distinct operation that draws no boundaries.

## Prerequisites

```bash
pip install squidpy scanpy scikit-image numpy pandas
# Image-based segmentation:
pip install cellpose          # v4 Cellpose-SAM generalist model
pip install stardist          # nuclear (star-convex) segmentation
# Transcript-based segmentation (CLI, on the molecule table):
#   Baysor (Julia binary) and proseg (Rust binary) are installed separately
```

## Quick Start

Tell your AI agent what you want to do:
- "Segment nuclei from my DAPI image with Cellpose"
- "Re-segment my Xenium data from the transcript table with Baysor"
- "Check whether my segmentation is over-segmenting or causing spillover"
- "Decide if my cell-by-gene matrix is trustworthy before cell typing"
- "Extract texture features under each Visium spot"

## Example Prompts

### Choosing a segmentation strategy
> "I have Xenium data with only DAPI -- should I use nucleus expansion, Cellpose, or Baysor, and what are the tradeoffs?"

> "My MERSCOPE run includes a membrane boundary stain. Which segmentation method makes use of it?"

### Running segmentation
> "Segment nuclei from this DAPI image and report how many cells were found"

> "Re-segment my CosMx data from the molecule table using proseg to recover immune cells"

### Judging trustworthiness
> "Some cells co-express EPCAM and CD3E -- is this a real hybrid state or segmentation spillover?"

> "Before I run ligand-receptor analysis, check whether short-range signals could be segmentation artifacts"

> "My transcripts-per-cell distribution is bimodal -- is this over- or under-segmentation?"

### Image features (not segmentation)
> "Extract summary and texture features under each Visium spot to help spatial-domain detection"

## What the Agent Will Do

1. Establish what boundary signal is available (DAPI only, membrane/boundary stain, or transcript table), because that gates the method choice.
2. Recommend a segmentation strategy from the tool ladder and explain its specific failure mode for the tissue at hand.
3. Run the chosen segmentation (Cellpose/StarDist on images, or Baysor/proseg on the molecule table) and build the cell-by-gene matrix.
4. QC the result: inspect transcripts-per-cell and cell-area distributions, flag impossible lineage co-expression, and run a contamination tool where available.
5. Report whether the derived matrix is provisional and which downstream analyses (rare states, hybrid states, short-range L-R) are most at risk.

## Tips

- The cell is a hypothesis, not an observation -- treat every cell-by-gene matrix from imaging data as provisional until QC'd.
- A membrane/boundary stain is the single highest-leverage change: it turns whole-cell boundary inference into measurement and beats swapping algorithms on DAPI-only data.
- Nucleus expansion assumes round, equal-sized cells; it is a baseline, not a solution, and is worst in dense or irregular tissue. Xenium cut its default expansion from 15 um to 5 um at XOA v2.0 -- an admission that 15 um over-assigned.
- Spillover is distance-dependent and strongest between adjacent heterotypic cells, so it fabricates spatially-structured false co-expression and the short-range co-localization that ligand-receptor tools detect -- a circularity to flag before trusting communication results.
- Standard doublet detectors (Scrublet, DoubletFinder) miss spatial doublets because those are neighbor merges, not the random pairs the detectors simulate.
- Segmentation-free methods (SSAM, ClusterMap) give the best cell-type maps but produce no cell objects -- per-cell composition, counts, and the neighbor graph are lost.
- Cellpose v4 (Cellpose-SAM) dropped the `channels=` argument and `diams` return value, and `diameter` is optional; older `model_type='nuclei'` code is v3.
- Squidpy `calculate_image_features` is NOT segmentation -- it summarizes the pixel patch under a spot and draws no boundaries.

## Related Skills

- spatial-transcriptomics/spatial-preprocessing - QC floors and non-gene-count normalization for the post-segmentation matrix
- spatial-transcriptomics/spatial-communication - ligand-receptor inference, where segmentation spillover fabricates short-range signal
- spatial-transcriptomics/spatial-proteomics - whole-cell segmentation on membrane markers for CODEX/IMC/MIBI
- imaging-mass-cytometry/cell-segmentation - segmentation for multiplexed-imaging proteomics
- spatial-transcriptomics/spatial-data-io - load the molecule table and the derived matrix

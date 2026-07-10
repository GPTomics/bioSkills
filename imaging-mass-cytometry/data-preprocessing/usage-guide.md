# IMC Data Preprocessing - Usage Guide

## Overview

Turns raw imaging mass cytometry (IMC) and MIBI acquisitions into compensated, variance-stabilized single-cell intensities. The fact that constrains every step: IMC pixels are integer ion COUNTS in a low-count Poisson regime, so every step (hot-pixel removal, spillover compensation, transformation, normalization) is a count-statistics problem where non-negativity is a hard physical constraint, spillover is spatial, and aggressive "cleaning" silently deletes real sparse biology.

## Prerequisites

```bash
# Python / CLI
pip install steinbock readimc numpy tifffile anndata

# R / Bioconductor (spillover compensation)
# BiocManager::install(c('CATALYST', 'cytomapper', 'imcRtools'))

# steinbock is also available as a Docker container
docker pull ghcr.io/bodenmillergroup/steinbock
```

## Quick Start

Tell your AI agent what you want to do:
- "Load my MCD files and extract per-channel TIFFs with hot-pixel removal"
- "Build the spillover matrix from my single-stain controls and compensate"
- "Compensate channel spillover at the pixel level before segmentation"
- "Arcsinh-transform my single-cell intensities with the right IMC cofactor"
- "Normalize across samples without breaking cross-sample comparability"

## Example Prompts

### Ingestion
> "Open this .mcd, list the acquisitions, and tell me which channels are metals versus antibody targets."

### Denoising
> "My CD20 and CD3 channels are very low-SNR. Should I run DeepSNiF on them, and which channels should I leave alone?"

### Spillover
> "I have single-stain control TXTs. Estimate the spillover matrix and compensate my images with NNLS. I care about spatial neighborhood analysis downstream, so compensate at the right level."

### Transformation
> "I came from suspension CyTOF and used cofactor 5. Is that right for IMC, and how should I normalize across 30 patient samples for differential abundance?"

## What the Agent Will Do

1. Ingest the multi-ROI MCD with readimc, carrying both metal (`channel_names`) and target (`channel_labels`) mappings.
2. Generate and respect a versioned `panel.csv` (the keep column filters and sorts the channel stack).
3. Remove hot pixels on raw counts with `--hpf` or IMC-Denoise DIMR; apply DeepSNiF only to low-SNR channels.
4. Estimate the spillover matrix from single-stain controls (per-event median) and compensate with NNLS -- at the pixel level before segmentation when spatial work is the endpoint, otherwise on cell means after segmentation.
5. Arcsinh-transform single-cell means with cofactor 1 and normalize against cohort-wide statistics, retaining raw counts.

## Tips

- The MCD is the source of truth; TXT is one ROI per file and is not emitted by Hyperion XTi.
- Never blanket-median IMC images -- it returns 0 for the isolated single-positive pixels that are often real signal.
- Compensation cannot rescue a saturated channel (>~5,000 dual counts); fix that at acquisition.
- If compensation produces negatives, the solver is wrong (use NNLS), not the concept; if it does nothing, the channel names are mismatched.
- State the arcsinh cofactor explicitly in any methods text -- there is no field-wide standard, so reproducibility depends on reporting it.
- Per-image percentile/min-max scaling is for pretty pictures, not quantitative comparison.

## Related Skills

- quality-metrics - reading the spillover matrix and gating channels before compensation
- cell-segmentation - segmentation runs on compensated nuclear/membrane channels
- phenotyping - consumes the arcsinh-transformed single-cell matrix
- flow-cytometry/compensation-transformation - suspension spillover and arcsinh background
- single-cell/preprocessing - AnnData conventions for the single-cell matrix

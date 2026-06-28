# Spatial Deconvolution - Usage Guide

## Overview

This skill estimates the cell-type composition of multi-cell spatial spots (Visium, Slide-seq, Stereo-seq bins, GeoMx regions) from an annotated scRNA-seq reference. It leads with the decision that gates everything: the resolution fork -- whether the platform produces mixtures that need deconvolving at all, or single cells that should be segmented and annotated instead. It then guides method choice (cell2location for absolute abundance vs RCTD/SPOTlight/stereoscope/SpatialDWLS for proportions), reference matching (the reference is the dominant determinant of the result), compositional handling of the sum-to-1 outputs, and reference-free STdeconvolve when no matched reference exists.

## Prerequisites

```bash
# cell2location (Python, pyro/scvi-tools; GPU recommended)
pip install cell2location

# Tangram (Python alternative)
pip install tangram-sc

# RCTD, SPOTlight, STdeconvolve (R)
# BiocManager::install(c('spacexr', 'SPOTlight', 'STdeconvolve'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Deconvolve my Visium spots into cell-type proportions using this scRNA-seq reference"
- "Estimate absolute cell numbers per spot with cell2location"
- "Is this Xenium data something I should deconvolve, or segment?"
- "I have no matched reference -- run a reference-free deconvolution"
- "Compare cell-type composition between tumor and normal regions"

## Example Prompts

### Deciding the regime
> "I have Xenium data -- should I run deconvolution on it?"

> "My platform is Slide-seqV2 with 10um beads -- deconvolve or reconstruct cells?"

### Running deconvolution
> "Deconvolve my Visium section using this annotated scRNA-seq reference and give me absolute cell abundances per spot"

> "Run RCTD in doublet-mode on my Slide-seq data"

### Reference and validation
> "A cell type I expect histologically is missing from my deconvolution output -- what happened?"

> "Check whether my estimated T-cell fraction tracks CD3D/CD3E expression across spots"

> "I don't have a matched scRNA-seq reference for this tissue -- deconvolve without one"

### Downstream
> "Compare cell-type composition between disease and control regions without the compositional artifact"

## What the Agent Will Do

1. Determine the resolution regime (deconvolve vs segment vs ambiguous) from the platform before choosing any tool.
2. Match the scRNA-seq reference to the tissue, condition, and ideally technology, and confirm it covers all expected cell types.
3. Choose a method by required output (absolute abundance vs proportions) and runtime, defaulting to cell2location or RCTD.
4. Train reference signatures, then map them onto the spatial data to estimate per-spot composition.
5. Transform compositional outputs with CLR/ILR before any downstream comparison.
6. Validate estimates against spatial marker-gene expression and, where useful, a reference-free STdeconvolve sanity check.

## Tips

- The resolution fork comes first. Running deconvolution on single-cell-resolution imaging data (Xenium, MERFISH, CosMx) is conceptually wrong -- it invents fractional mixtures inside cells that are already pure. Segment and annotate those instead.
- The reference IS the result. A cell type missing from the reference is silently reassigned to its nearest present neighbor, with no error flag -- confident, wrong, and invisible.
- Match the reference to tissue AND condition. A healthy reference mis-estimates activated-immune and malignant states whose expression has shifted.
- Reference quality beats method choice. Benchmarks show a plain NNLS outperforms half the dedicated methods, and reference quality swings results more than the algorithm.
- Rare-type fractions (below a few percent) are the least reliable numbers -- corroborate with spatial markers before believing them.
- cell2location alone returns absolute cell abundance; the others return proportions only. Pick by the biological question at hand.
- Outputs sum to 1. Use CLR/ILR (or ALDEx2/scCODA) for downstream comparison, never naive per-type t-tests.
- Two NB-regression methods agreeing is pseudo-replication, not validation. Perturb the reference and use an orthogonal modality to build confidence.

## Related Skills

- spatial-domains - group spots into tissue regions; a domain is a region, not a cell type
- high-resolution-binning - the ambiguous near-single-cell regime (Visium HD, Stereo-seq, Slide-seq)
- image-analysis - segment cells from imaging platforms, where deconvolution does not apply
- single-cell/cell-annotation - annotate imaging cells and label the scRNA-seq reference
- single-cell/preprocessing - build a clean reference before deconvolving
- spatial-preprocessing - QC and normalize spatial data before deconvolution
- spatial-visualization - map proportions and abundances onto the tissue

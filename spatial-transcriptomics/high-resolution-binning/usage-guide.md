# High-Resolution Binning

## Overview

High-resolution spatial platforms capture transcripts at or below the size of a single cell: Visium HD records 2um square bins, Stereo-seq records ~220nm DNB spots, and Slide-seqV2 records 10um beads. Because each capture unit is smaller than (or roughly the size of) one cell, the analytical task is the INVERSE of deconvolution -- aggregate sub-cellular fragments UP into whole cells rather than mixing a multi-cell spot DOWN into cell-type fractions. This skill covers the central bin-size dilemma (2um bins are too sparse to cluster, but binning to 8/16um re-creates the multi-cell mixture), the decision between morphology-driven cell reconstruction (Bin2cell -- segment nuclei on a registered image, then assign 2um bins to nuclei) and fixed-bin aggregation, and how that decision is set by the platform and whether a co-registered cell image exists.

## Prerequisites

- Python with scanpy, anndata, and numpy
- Bin2cell for Visium HD morphology-driven reconstruction: `pip install bin2cell`
- StarDist or Cellpose for nucleus segmentation during reconstruction: `pip install stardist` (needs TensorFlow) or `pip install cellpose`
- spatialdata and spatialdata-io for loading high-resolution outputs: `pip install spatialdata spatialdata-io`
- squidpy for downstream spatial analysis: `pip install squidpy`
- For Slide-seqV2 bead deconvolution as an alternative: RCTD (R package spacexr); see spatial-deconvolution

```bash
pip install bin2cell scanpy anndata spatialdata spatialdata-io squidpy
pip install stardist   # or: pip install cellpose
```

## Quick Start

Tell your AI agent what you want to do:
- "Turn my Visium HD 2um bins into single cells"
- "Reconstruct cells from my Visium HD H&E image with Bin2cell"
- "Should I analyze my Visium HD data at 2um, 8um, or 16um?"
- "Aggregate my Slide-seqV2 beads into cell-scale profiles"
- "My Stereo-seq bins span multiple cells -- how do I get to single cells?"
- "Is this a binning problem or a deconvolution problem?"

## Example Prompts

### Choosing the approach
> "I have Visium HD output with a registered H&E image. Walk me through reconstructing single cells instead of just using the 8um bins, and explain why the 8um bins are not cells."

> "My Slide-seqV2 puck has no cell image. What is the right way to get cell-scale profiles, and when should I deconvolve the beads instead?"

### Bin size and the mixture dilemma
> "Explain the trade-off between 2um, 8um, and 16um Visium HD bins for my analysis, and which one I should cluster on."

> "I tried clustering my 2um bins and got noise. What went wrong and what should I do instead?"

### Recognizing the inverse-of-deconvolution framing
> "Someone told me to deconvolve my Visium HD 2um bins. Is that correct?"

> "Is reconstructing cells from Stereo-seq bins the same problem as deconvolving Visium spots?"

## What the Agent Will Do

1. Identify the platform and the native capture-unit size, and place the dataset in the AMBIGUOUS regime of the resolution fork (unit smaller than or comparable to one cell).
2. Check whether a cell-resolution image is co-registered to the capture coordinates, the single fact that decides reconstruction vs aggregation.
3. For Visium HD with a registered image, run morphology-driven reconstruction: scale the image, destripe the bins, segment nuclei with StarDist or Cellpose, assign 2um bins to nuclei, expand to capture cytoplasmic bins, and collapse bins per nucleus into cell-level profiles.
4. For Slide-seqV2 or imageless data, aggregate to a cell-scale grid (chosen in microns near one cell diameter) or route the beads to bead-level deconvolution when they clearly mix cell types.
5. Warn against the coarse-binning trap (8um bins still span ~2 cells) and against deconvolving sub-cellular fragments (which invents mixtures that do not exist).
6. Apply cell-scale QC to the reconstructed cells (filter on bins-per-cell, the way single-cell QC filters on UMIs) before handing off to clustering and annotation.

## Tips

- The capture unit is a fragment of one cell, not a mixture of several -- aggregate UP, never deconvolve DOWN. Deconvolution is for spots LARGER than a cell.
- 2um bins are too sparse to cluster (single-digit median UMIs); do not feed raw bins to PCA or Leiden.
- Binning to 8um or 16um does not solve sparsity -- it re-creates the multi-cell mixture, putting the data back in the deconvolution regime with a reference required.
- The discriminating factor is the registered cell image, not the platform name: Visium HD without a usable image behaves like Slide-seqV2; Stereo-seq with a clean registered stain behaves like Visium HD.
- Verify the morphology image is registered to the bin coordinate frame before trusting reconstruction; a misregistered image assigns bins to the wrong nuclei with no error.
- Filter reconstructed cells on bins-per-cell; a cell built from very few bins is a low-confidence profile.
- Run destripe before segmentation on Visium HD -- per-row/per-column total-count striping otherwise biases both nuclei and counts.
- Bins assigned to no nucleus are inter-cellular space; dropping them is correct, forcing them into a cell fabricates expression.

## Related Skills

- spatial-deconvolution - the resolution fork that routes the AMBIGUOUS regime here; deconvolution is the opposite (mix DOWN) geometry
- image-analysis - nucleus and cell segmentation (StarDist, Cellpose) underpinning morphology-driven reconstruction
- spatial-data-io - load Visium HD PARQUET bin positions and the registered image
- spatial-preprocessing - QC and normalize the reconstructed cells at cell scale, not bin scale
- single-cell/cell-annotation - annotate the reconstructed cells with markers or label transfer
- single-cell/clustering - cluster reconstructed cells once they carry cell-scale signal

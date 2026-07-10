# Spatial Data I/O - Usage Guide

## Overview

This skill loads spatial transcriptomics data from imaging/in-situ platforms (Xenium, MERSCOPE/MERFISH, CosMx) and sequencing/capture platforms (Visium, Visium HD, Slide-seq/Curio, Stereo-seq) into AnnData or SpatialData objects. It centers on the I/O decision that governs everything downstream: imaging platforms emit two distinct primary objects -- a per-transcript molecule table (the re-segmentable source of truth) and a segmentation-derived per-cell matrix (quality-filtered, inherits all segmentation error) -- while capture platforms emit only a spot/bin matrix with no molecule table. It also covers picking the right reader per platform and keeping coordinate frames (pixel vs micron) registered to histology.

## Prerequisites

```bash
pip install spatialdata spatialdata-io squidpy scanpy anndata
```

## Quick Start

Tell your AI agent what you want to do:

- "Load my Visium data from the Space Ranger output"
- "Read my Xenium experiment and keep the transcript molecule table"
- "Load my MERSCOPE / MERFISH data"
- "Import my CosMx flat files"
- "Load my Visium HD bins"
- "Which reader do I use for Slide-seq?"

## Example Prompts

### Sequencing / capture platforms
> "Load my 10X Visium Space Ranger output and tell me whether the coordinates are in pixels or microns."

> "Read my Visium HD bundle at the 8 micron bin and remind me whether bins are single cells."

> "Load my Slide-seq / Curio beads; do they have a transcript molecule table?"

### Imaging / in-situ platforms
> "Load my Xenium experiment as a SpatialData object and give me both the molecule table and the cell matrix."

> "Read my MERSCOPE output -- I want the per-transcript detections, not just the cell-by-gene matrix."

> "Import my CosMx data and check whether a transcript file is present."

### Reasoning about objects
> "Is the Xenium cell-by-gene matrix the same quality filter as the transcript table?"

> "I have Visium spots -- where is the molecule table?"

## What the Agent Will Do

1. Determine the platform class (imaging/in-situ vs sequencing/capture), which decides what objects exist and whether a molecule table is present.
2. Select the matching reader (`spatialdata_io.{xenium, merscope, cosmx, visium, visium_hd, curio, stereoseq}` or `squidpy.read.{visium, vizgen, nanostring}`), avoiding the deprecated `scanpy.read_visium`.
3. Load into SpatialData (imaging, to retain the molecule table, shapes, and images) or AnnData (spot matrices), and report element names.
4. Distinguish the per-transcript molecule table from the segmentation-derived, quality-filtered cell matrix, flagging the matrix as provisional.
5. Confirm the coordinate frame and units before any graph-building or histology overlay.

## Tips

- The imaging cell-by-gene matrix is a DERIVED product of segmentation and is usually quality-filtered (Xenium keeps Q>=20); the molecule table keeps everything and is the only object that lets the data be re-segmented. Do not trust the matrix as ground truth.
- Spot platforms (Visium, Slide-seq, Stereo-seq) have NO molecule table -- a spot/bead/bin is mini-bulk over several cells. Do not go looking for transcript detections that do not exist.
- `squidpy.read` provides only `visium`, `vizgen`, and `nanostring`. For Xenium, Slide-seq/Curio, Stereo-seq, and Visium HD, use `spatialdata_io`.
- There is no `merfish` reader -- `spatialdata_io.merscope` handles MERFISH and MERSCOPE both.
- Visium HD tissue positions are PARQUET, not CSV; the Visium `tissue_positions.csv` gained a header at Space Ranger v2.0 (readers handle both).
- `obsm['spatial']` units are not universal: Visium stores full-res pixels (scale with `tissue_hires_scalef`), most imaging readers place a micron "global" frame. Check before building a neighbor graph or overlaying transcripts.
- A SpatialData `.zarr` store is a DIRECTORY -- write it to a scratch path and `rm -rf` it; never commit it.

## Related Skills

- spatial-preprocessing - QC floors and normalization that differ by platform class after loading
- image-analysis - re-segment the molecule table; the cell matrix is a segmentation hypothesis
- spatial-deconvolution - recover cell-type proportions from spot mixtures that have no molecule table
- high-resolution-binning - bin/segment-up sub-cellular Visium HD and Stereo-seq captures
- spatial-visualization - plot spots vs imaging FOVs with the correct coordinate frame
- single-cell/data-io - non-spatial scRNA-seq loading for the deconvolution/label-transfer reference

# Cell Segmentation - Usage Guide

## Overview

Delineates single cells in multiplexed IMC/MIBI images so that averaging channels inside each mask yields a single-cell expression table. What dominates everything downstream: segmentation is the single largest, irreversible source of error in the experiment -- under-segmentation fabricates impossible double-positive cell types while over-segmentation quietly miscounts cells and corrupts spatial statistics, and both are confounded with lateral signal spillover across cell boundaries.

## Prerequisites

```bash
# steinbock orchestrates DeepCell/Mesmer and Cellpose
pip install steinbock deepcell cellpose scikit-image numpy tifffile

# Docker container alternative
docker pull ghcr.io/bodenmillergroup/steinbock
```

## Quick Start

Tell your AI agent what you want to do:
- "Segment whole cells with Mesmer using my DNA and membrane channels"
- "Build a membrane channel that covers all my cell types"
- "Fall back to nuclear segmentation plus expansion -- my membrane staining is poor"
- "Check my segmentation for impossible double-positive cells"
- "Run steinbock cellpose segmentation and extract per-cell mean intensities"

## Example Prompts

### Whole-cell segmentation
> "Run Mesmer whole-cell segmentation on this IMC image. My resolution is 1 um/pixel -- set image_mpp correctly."

### Membrane-channel decision
> "Which markers should I sum into the membrane channel so I don't under-segment my B cells and macrophages?"

### Fallback
> "My only good membrane marker is CD45 and it only stains immune cells. What is the least-biased way to get whole-cell masks?"

### Evaluation
> "I have a CD3+CD20+ cluster. Is it real, and how do I tell if it is a segmentation artifact?"

## What the Agent Will Do

1. Choose a segmenter from the decision tree (Mesmer first for whole-cell; nuclei+expansion when membrane staining is weak; ilastik+CellProfiler for transparency).
2. Build the nuclear and summed-membrane inputs, covering all cell types in the membrane sum.
3. Run segmentation with the true `image_mpp`/`diameter`, then aggregate per-cell MEAN intensities.
4. Apply lateral-spillover compensation (REDSEA) after segmentation when membrane double-positives are heavy.
5. Evaluate on downstream proxies -- impossible-co-expression rate, count/density sanity, positive-fraction stability across segmenters -- not F1/IoU alone.

## Tips

- A biologically-impossible co-expressing population is a segmentation diagnosis until proven otherwise.
- Pass the TRUE acquisition resolution to Mesmer's `image_mpp` (~1.0 IMC); the default rescales cells to the wrong size.
- Set Cellpose `diameter` deliberately at IMC scale, or use `cpsam` which removes the diameter dependence.
- Use exclusive (watershed/`expand_labels`) expansion, never free dilation -- free dilation double-counts boundary pixels.
- Channel compensation goes before aggregation (on pixels); lateral-spillover REDSEA goes after segmentation (on cells).
- F1/IoU on sparse annotations is optimistic and off-target; audit dense regions, not easy ones.

## Related Skills

- data-preprocessing - channel spillover compensation precedes segmentation
- phenotyping - consumes the single-cell mask and intensities; double-positives diagnose segmentation
- spatial-analysis - over-segmentation corrupts neighborhood statistics
- quality-metrics - segmentation QC metrics and the impossible-co-expression monitor
- interactive-annotation - overlay masks on channels to audit boundaries

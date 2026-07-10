# Interactive Annotation - Usage Guide

## Overview

Labels cells and QCs segmentation by looking at the image with masks overlaid. Why this step matters: annotation is the bridge between the raw pixel image and the single-cell table, and it is the only QC step that can catch the segmentation, spillover, and morphology artifacts that summary statistics on the cell table structurally cannot show. Manual labels are not ground truth (expert-vs-expert concordance is only ~86%), and the display contrast limit is a positivity threshold in disguise.

## Prerequisites

```bash
pip install napari napari-imc scikit-learn numpy
# GUI backend
pip install pyqt5

# R viewers (optional): BiocManager::install(c('cytomapper', 'cytoviewer'))
# Mantis Viewer is a separate Electron app (CANDELbio / Parker Institute)
```

## Quick Start

Tell your AI agent what you want to do:
- "Open my .mcd in napari with the segmentation mask overlaid"
- "Paint my clusters back onto the tissue and tell me which are artifacts"
- "Build a class-balanced training set across my patients"
- "Set up per-cell type annotation with fixed contrast limits"
- "Check whether my CD3+CD20+ cells are merged segments on the image"

## Example Prompts

### Setup
> "Load this .mcd and the cell mask in napari, overlay outlines on DNA and a summed membrane channel, and fix the contrast so positivity is consistent."

### Cluster validation
> "Color my masks by Leiden cluster and tell me which clusters scatter as salt-and-pepper haze along boundaries -- those are the artifact clusters."

### Ground truth
> "I need to train CellSighter. How many cells per class should I label, how do I handle my rare tumor population, and across how many patients?"

### Artifact check
> "Overlay the mask and the CD3 and CD20 channels on my suspected double-positive cells and tell me if the two signals come from two abutting nuclei."

## What the Agent Will Do

1. Open the raw `.mcd` with napari-imc and overlay mask outlines on nuclear + summed-membrane channels at fixed contrast limits.
2. Paint data-driven clusters back onto tissue and flag spatially-incoherent clusters as artifacts.
3. Set up per-cell annotation (points layer with categorical labels) or mask QC (labels layer).
4. Build a class-balanced ground-truth set, over-sampling rare types across multiple patients.
5. Route high-uncertainty cells (Astir probability, CellSighter score) back for human review.

## Tips

- Distrust any cluster, double-positive, or rare type until it has been seen on the image with its mask.
- Core napari cannot open `.mcd`; install the napari-imc plugin.
- Fix and record contrast limits -- auto-scaling per image silently moves the positivity threshold.
- Annotate on tissue with neighbors visible, not from a biaxial/UMAP scatter that strips context.
- Do not chase classifier accuracy above the ~86% inter-annotator ceiling; use consensus labels for evaluation.
- Channel-spillover compensation does not handle lateral (neighbor) spillover -- catch that visually.

## Related Skills

- cell-segmentation - mask overlay is the irreplaceable segmentation QC
- phenotyping - annotation supplies labels/priors and confirms clusters are real
- quality-metrics - the image catches artifacts that table statistics cannot
- data-preprocessing - contrast/transform choices mirror preprocessing thresholds
- spatial-analysis - spatial coherence of a painted cluster validates it

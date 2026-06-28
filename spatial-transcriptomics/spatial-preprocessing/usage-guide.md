# Spatial Preprocessing - Usage Guide

## Overview

This skill covers quality control, filtering, and normalization for spatial transcriptomics data, with the decisions forked by platform class. Sequencing/spot data (Visium, Visium HD, Slide-seq, Stereo-seq) and imaging/in-situ data (Xenium, MERSCOPE/MERFISH, CosMx, seqFISH) need different QC floors and, critically, different normalization. The central reframe is that in spatial data library size carries biology rather than pure technical depth, so the scRNA reflex of dividing it out removes real spatially-structured signal and degrades downstream domain detection. The skill gives per-platform QC floors, the negative-control false-discovery metric for imaging specificity, and a normalization decision table that steers skewed targeted panels toward cell-area normalization instead of Pearson residuals.

## Prerequisites

```bash
pip install squidpy scanpy anndata spatialdata
```

## Quick Start

Tell your AI agent what you want to do:
- "QC and normalize my Visium data"
- "Set QC floors for my Xenium cells without deleting real cells"
- "Compute the negative-control false-discovery rate for my MERFISH panel"
- "Normalize my imaging data by cell area instead of library size"
- "Show QC metrics on the tissue to check for spatial artifacts"

## Example Prompts

### QC floors by platform
> "My Xenium cells have about 100 transcripts each -- what count floor should I use so I do not delete real lymphocytes?"

> "Run QC on my Visium section using tissue-appropriate UMI, gene, and mito thresholds and show me where the cut spots fall."

### Specificity and controls
> "Compute the negative-control-probe and blank-barcode false-discovery rate for my imaging panel and drop the controls before clustering."

### Normalization decision
> "Should I normalize my spatial data at all, given that library size tracks cellularity? Pick a method for my skewed targeted panel."

> "Normalize my Xenium data by segmented cell area rather than total counts and explain why Pearson residuals would not fix the panel-skew bias."

### Spatial QC inspection
> "Map total counts and genes per cell onto the tissue to check for an edge or permeabilization gradient before I threshold."

## What the Agent Will Do

1. Identify the platform class (spot vs imaging) to fork every threshold decision.
2. Annotate negative controls and mito genes, then compute QC metrics.
3. Compute the negative-control false-discovery rate as the imaging specificity metric.
4. Map QC metrics onto tissue coordinates to catch spatial artifacts before cutting.
5. Apply platform-appropriate floors that do not preferentially delete small real cells.
6. Choose a normalization method from the decision table -- cell-area/volume for skewed imaging panels, deliberate log1p or Pearson residuals for whole-transcriptome spot data.

## Tips

- An scRNA `min_counts=500` floor deletes nearly every imaging cell (their vectors are tens-to-low-hundreds of transcripts) -- use ~10 transcripts/cell, a community convention, not a vendor spec.
- "High genes per cell = doublet" is meaningless for imaging: genes/cell is ceilinged at the panel size.
- Mito-% QC is usually impossible on imaging panels because mito genes are off-panel -- do not block the pipeline waiting for it.
- Library size carries biology in spatial data: on Visium it tracks cells-per-spot and cellularity, on imaging it tracks cell size/area -- dividing it out blurs spatial domains.
- Pearson residuals and SCTransform do NOT rescue a skewed imaging panel; they are still gene-count-based. Prefer cell-area/volume normalization or SpaNorm.
- Negative-control FDR (mean-control / mean-real-gene) is the only native false-discovery proxy for imaging -- compute it, then drop controls before clustering.
- Always inspect QC spatially: a smooth gradient across the section is a technical artifact, not biology, and a violin plot hides it.

## Related Skills

- spatial-data-io - load Visium/Xenium/MERFISH and reach the molecule table vs the segmentation-derived matrix
- image-analysis - segment cells from imaging data, the upstream error source that sets imaging QC and cell area
- spatial-deconvolution - the next step for spot data, where library size and reference choice decide proportions
- single-cell/preprocessing - the scRNA QC/normalization baseline these thresholds deliberately depart from

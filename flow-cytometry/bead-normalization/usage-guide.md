# Bead Normalization - Usage Guide

## Overview
Bead-based normalization corrects CyTOF instrument sensitivity drift using EQ four-element beads as a physical internal standard, while cross-batch normalization (CytoNorm) harmonizes staining/acquisition batches using a shared reference sample. The key expert distinction this skill encodes: these are two different layers - beads fix within/between-run drift, CytoNorm fixes batch effects - and they are not interchangeable. It also covers why the anchor/reference sample is the design element everything depends on, why normalization is per-cluster with many quantiles, and the over-correction risk that argues for modeling batch in the statistical design rather than cleaning it out of the data.

## Prerequisites
```bash
# R/Bioconductor
BiocManager::install(c('CATALYST', 'flowCore'))
# CytoNorm (cross-batch)
# remotes::install_github('saeyslab/CytoNorm')
```

## Quick Start
Tell your AI agent what you want to do:
- "Normalize my CyTOF run with EQ beads and remove the bead events"
- "Harmonize my multi-batch study using my reference PBMC sample"
- "Plot bead signal over acquisition time to check for drift"
- "Should I normalize batch out or model it in my statistics?"

## Example Prompts
### EQ-bead drift correction
> "Run CATALYST normCytof on this CyTOF SCE using the DVS EQ bead masses, correct the drift, remove the beads, and tell me what fraction were bead events."
> "Plot the Ce140 bead-channel median over Time before and after normalization to confirm the drift is gone."

### Cross-batch
> "Train CytoNorm on my reference sample that was run in all four batches, then normalize the experimental samples - and run testCV first so I know the per-cluster splines are safe."
> "One of my batches doesn't have the reference aliquot - what are my options?"

### Strategy
> "I'm going to do differential abundance downstream - should I CytoNorm the data or put batch in the diffcyt design? Explain the over-correction tradeoff."

## What the Agent Will Do
1. Apply EQ-bead normalization first (on raw counts) to correct sensitivity drift, returning the cleaned SCE from `normCytof()$data`.
2. For multi-batch data, confirm an anchor/reference sample exists in every batch.
3. Run `testCV()` to verify FlowSOM clustering is batch-stable, then train and apply CytoNorm per-cluster.
4. Visualize bead/marker signal before and after to confirm correction without over-correction.
5. Recommend modeling batch in the diffcyt design for inference, reserving normalization for display/clustering.

## Tips
- Bead normalization (drift) and CytoNorm (batch) are different layers - do bead first, on raw counts; CytoNorm after, on transformed scale.
- `normCytof()` returns a list - the normalized SCE is in `$data`; `beads="dvs"` is the EQ four-element mass set.
- CytoNorm requires the same anchor/reference sample in every batch; without it, that batch can't be safely normalized.
- Run `testCV()` first - if per-cluster CV is high the FlowSOM model is batch-unstable and splines will distort (fall back to `nClus=1`).
- Per-cluster + ~99 quantiles handles cell-type-specific, non-linear drift; a single global shift reintroduces distortion.
- Over-normalization can erase real biology; for inference prefer modeling batch in the design, normalizing only for visualization.
- Fully confounded batch (batch == condition) cannot be rescued by any normalization.

## Related Skills
- cytometry-qc - EQ-bead-median-vs-Time is the primary CyTOF drift readout
- doublet-detection - Remove doublets before normalization
- compensation-transformation - The transform scale CytoNorm operates on
- clustering-phenotyping - Cluster across normalized batches
- differential-analysis - Model batch in the design rather than over-cleaning
- experimental-design/batch-design - Anchor/reference-sample design and batch strategy

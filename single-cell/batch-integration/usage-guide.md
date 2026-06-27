# Batch Integration - Usage Guide

## Overview

Batch integration learns a shared low-dimensional representation that mixes technical batches while preserving biological cell states. The central tension is that batch-mixing and biological-signal preservation cannot be jointly maximized: over-correction silently erases rare cell types and continuous gradients, and under confounded designs batch and biology are unidentifiable. This skill resolves which method to use (Harmony, scVI/scANVI, Seurat CCA/RPCA, fastMNN, Scanorama, BBKNN), how strongly to correct, when not to integrate, how to score integration with scIB metrics without gaming them, and why corrected expression must never feed differential expression.

## Prerequisites

**Python:**
```bash
pip install scanpy harmonypy scvi-tools scanorama bbknn scib-metrics
```

**R:**
```r
install.packages(c('Seurat', 'harmony'))
BiocManager::install('batchelor')
```

## Quick Start

Tell the agent what to do:
- "Integrate my samples and remove batch effects"
- "Pick an integration method for a large atlas with many batches"
- "Is my batch confounded with condition - should I integrate at all?"
- "Score my integration and check I did not over-correct"

## Example Prompts

### Integration
> "Merge my samples and run Harmony, then cluster on the corrected embedding"

> "Use scVI to integrate a large multi-batch atlas from raw counts"

> "Run Seurat v5 RPCA integration because the datasets only partly overlap"

### Method Choice
> "I have a few same-protocol samples - which method and should I even integrate?"

> "Some cells are labeled - would scANVI preserve biology better here?"

### Diagnosis
> "A rare population disappeared after integration - is this over-correction?"

> "My batch tracks my treatment - can integration separate them?"

### Scoring
> "Compute kBET, iLISI, and cell-type silhouette and tell me if biology was kept"

> "Benchmark Harmony vs scVI vs scANVI with scib-metrics and pick the most robust"

## What the Agent Will Do

1. Visualize uncorrected data first and cross-tabulate clusters x batch x condition to rule out confounding.
2. Decide whether to integrate at all (skip for confounded designs, well-mixed technical replicates, per-sample analyses).
3. Pick a method by dataset size and design (Harmony/Seurat for simple/small, scVI/scANVI for large/complex/label-rich, Scanorama/RPCA for partial overlap, BBKNN for speed).
4. Run integration on a chosen strength, keeping the uncorrected embedding for comparison.
5. Score batch-mixing and bio-conservation separately (scIB composite = 0.6 bio + 0.4 batch), checking rare populations before/after.
6. Use the integrated embedding for clustering/visualization only; run DE on uncorrected counts via pseudobulk.

## Tips

- **Visualize before correcting** - integration is a bias-variance trade; do not integrate on reflex.
- **Confounded design is unfixable computationally** - if each condition is its own batch, no method separates them; redesign with multiplexing.
- **Over-correction is the silent failure** - it erases rare types and gradients while batch metrics improve; the cells worth finding are the ones most at risk.
- **No universal winner** - Harmony/Seurat for simple/small, scVI/scANVI for large/complex/label-rich (scIB, Luecken 2022); run 2-3 candidates and score.
- **Batch metrics are gamed by over-correction** - kBET/iLISI are maximized by destroying structure; always pair with bio metrics (ASW-celltype, cLISI, ARI).
- **scVI latent is not "biology minus batch"** - it is an entangled coordinate system; do not interpret individual dimensions.
- **Never run DE on corrected expression** - use uncorrected log-normalized counts; pseudobulk per sample x cell type for cross-condition DE.
- **Reference mapping is closed-world** - projected labels look confident even for novel states; inspect mapping uncertainty.

## Related Skills

- preprocessing - QC and normalization that must precede integration
- clustering - Cluster on the integrated embedding, not on raw PCA
- cell-annotation - Reference mapping and label transfer after integration
- single-cell/multimodal-integration - Joint analysis across modalities (distinct from batch integration)
- single-cell/differential-abundance - Test whether composition shifts across conditions after integration
- differential-expression/deseq2-basics - Pseudobulk DE on uncorrected counts per cell type
- data-visualization/dimensionality-reduction-plots - Before/after UMAP comparison figures

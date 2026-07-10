# Reference Mapping and Label Transfer Usage Guide

## Overview

Annotate a query single-cell dataset by projecting it onto a pre-trained reference atlas and transferring cell-type labels, using scArches surgery (scVI/scANVI), Symphony, Azimuth, CellTypist, scPoli, popV, or fine-tuned foundation models. The central discipline this skill enforces is that mapping is a projection, not an annotation oracle: every transferred label must be gated on an out-of-distribution / label-transfer-uncertainty signal, which is a different quantity from the classifier probability.

## Prerequisites

```bash
pip install scvi-tools scanpy anndata scikit-learn celltypist
```

Conceptual prerequisites: the query must be QC'd (doublets removed, ambient corrected, MT/count filters) to the reference's standard before mapping, and its genes must align to the reference feature set.

## Quick Start

Tell your AI agent what you want to do:
- "Map my single-cell data to the Human Lung Cell Atlas and flag cells that do not belong"
- "Transfer cell-type labels from my reference to my query and gate them on uncertainty"
- "Pick the right mapping method for a clinical, reproducible annotation pipeline"
- "My query likely has disease states the reference lacks; how do I avoid confident wrong labels?"

## Example Prompts

### Method choice

> "I need a deterministic, CPU-only, reproducible annotation pipeline for a clinical setting. Which mapping method should I use and why?"

> "Should I use scGPT zero-shot embeddings to annotate my query, or scANVI surgery?"

### Label transfer with uncertainty

> "Map my query to the reference scANVI model with scArches surgery, then set low-confidence cells to Unknown using weighted-kNN transfer uncertainty, not the softmax."

> "Which cells in my query are out-of-distribution relative to the reference? They may be novel cell types or disease states."

### Trajectories and novel biology

> "My query is a differentiation trajectory. Will scANVI collapse it onto discrete reference labels, and should I use unsupervised scVI surgery instead?"

### Evaluation

> "I have a high scIB integration score. Does that mean my transferred labels are correct?"

## What the Agent Will Do

1. Choose a mapping method from the taxonomy based on reproducibility needs, novel-biology risk, and whether an embedding is required
2. QC the query to the reference standard and align genes with `prepare_query_anndata`
3. Run surgery (or linear/classifier mapping) with frozen reference weights and `weight_decay=0.0`
4. Transfer labels and the latent embedding
5. Compute a label-transfer uncertainty / OOD signal and gate low-confidence cells to "Unknown"
6. Evaluate with held-out labels, marker sanity checks, and OOD on spiked-in unseen types -- not integration score alone

## Tips

- The softmax max from `predict(soft=True)` answers "which label," not "does it belong" -- gate on weighted-kNN uncertainty (HLCA threshold 0.2) or a Mahalanobis/ensemble OOD signal
- `weight_decay=0.0` and frozen weights are not optional: they keep the shared latent fixed so queries mapped to the same atlas stay comparable
- `prepare_query_anndata` is mandatory; skipping it silently corrupts the encoder input rather than raising an error
- For trajectories or suspected novel states, prefer unsupervised scVI surgery and annotate the query independently
- CellTypist expects log1p of CP10k input and is a classifier, not a mapper (no shared embedding, no batch handling)
- Zero-shot foundation-model embeddings are not a default; they underperform scVI/Harmony/HVG-PCA on integration (Kedzierska 2025)
- cellxgene-census provides direct access to CZI atlases (not `scvi.data.cellxgene()`)

## Related Skills

- single-cell/preprocessing - QC, normalization, and HVG selection before mapping
- single-cell/markers-annotation - Manual marker-based cluster annotation without a reference
- single-cell/batch-integration - Integrating datasets without a labeled reference
- single-cell/doublet-detection - Removing doublets that map to spurious intermediates
- differential-expression/de-results - Pseudobulk validation of mapping-derived populations

# Transfer Learning Usage Guide

## Overview

Map query single-cell datasets to pre-trained reference atlases using scArches with scVI and scANVI models, enabling rapid cell type annotation without retraining on combined data.

## Prerequisites

```bash
pip install scvi-tools scanpy anndata
```

## Quick Start

Tell your AI agent what you want to do:
- "Map my single-cell data to the Human Lung Cell Atlas"
- "Transfer cell type labels from my reference to query dataset"
- "Use scANVI to predict cell types in my new dataset"
- "Get prediction confidence scores for transferred labels"

## Example Prompts

### Reference Mapping

> "I have a pre-trained scVI model from my reference atlas. Map my new query dataset to the same latent space using scArches surgery."

> "Transfer my query single-cell data to the reference embedding without retraining the full model."

### Label Transfer

> "Use scANVI to transfer cell type labels from my annotated reference to my unlabeled query dataset. Show prediction confidence."

> "Which cells in my query have low-confidence label predictions? They might be novel cell types."

### Visualization

> "Create a UMAP showing my query cells embedded with the reference atlas. Color by dataset and cell type."

## What the Agent Will Do

1. Load pre-trained reference model and query data
2. Subset query to reference genes
3. Set up query AnnData with reference schema
4. Run surgical fine-tuning (weight_decay=0.0)
5. Transfer latent representation and/or labels
6. Report prediction confidence

## Tips

- Query must have same genes as reference (subset before mapping)
- Use `weight_decay=0.0` for surgical training to preserve reference structure
- scANVI transfers labels; scVI only transfers embedding
- Low-confidence predictions may indicate novel cell types or poor data quality
- For large atlases, download pre-trained models from CellxGene or similar
- cellxgene-census provides direct access to CZI atlases (not scvi.data.cellxgene())

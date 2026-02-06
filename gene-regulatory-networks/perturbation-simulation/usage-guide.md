# Perturbation Simulation - Usage Guide

## Overview

Simulate transcription factor perturbation effects on cell state using CellOracle. Constructs a GRN by combining a base GRN from chromatin accessibility data with scRNA-seq expression, then predicts how cell identities shift when TFs are knocked out or overexpressed. Useful for prioritizing TF perturbation experiments and identifying driver TFs for cell fate transitions.

## Prerequisites

```bash
pip install celloracle scanpy matplotlib
```

CellOracle requires:
- scRNA-seq data (preprocessed AnnData with clustering and UMAP)
- Base GRN from accessible chromatin regions (scATAC, bulk ATAC, or published data)

Note: CellOracle does NOT require paired multiome data. Any source of accessible regions can provide the base GRN.

## Quick Start

Tell your AI agent what you want to do:
- "Simulate a TF knockout and show how cell states change"
- "Which transcription factors drive differentiation in my data?"
- "Predict what happens if I knock out GATA1 in my hematopoietic cells"
- "Screen multiple TFs to find drivers of cell fate"

## Example Prompts

### Single TF Perturbation
> "I have scRNA-seq data and ATAC peaks. Simulate knocking out PAX5 and show the predicted cell state shifts."

> "Overexpress CEBPA in silico and visualize which cells are most affected."

### TF Screening
> "Screen these 10 TFs for their impact on cell fate: GATA1, SPI1, CEBPA, PAX5, TCF7, RUNX1, EBF1, IRF4, FOXP3, TBX21."

> "Rank TFs by their predicted effect on my progenitor cell cluster."

### Visualization
> "Create a quiver plot showing the direction of cell state change after GATA1 knockout."

> "Show which cells are most affected by this perturbation on the UMAP."

## What the Agent Will Do

1. Scan accessible regions for TF binding motifs to build the base GRN
2. Import scRNA-seq data and base GRN into CellOracle
3. Fit GRN models per cell type using regularized regression
4. Simulate TF perturbation (knockout or overexpression)
5. Calculate predicted cell state shifts on the embedding
6. Visualize results as quiver plots and gradient plots

## Tips

- **Base GRN source** - Can come from scATAC-seq, bulk ATAC-seq, or published chromatin data; does not require paired multiome
- **Pre-built base GRNs** - CellOracle provides pre-computed base GRNs for common mouse and human tissues
- **n_propagation** - Controls how far perturbation effects propagate through the GRN; default of 3 is usually sufficient
- **Systematic screening** - Loop over candidate TFs to rank them by predicted impact; focus on TFs expressed in your cell types
- **Validation** - Compare predictions with actual Perturb-seq data or published KO phenotypes when available

## Related Skills

- scenic-regulons - TF regulon inference from scRNA-seq with pySCENIC
- multiomics-grn - Enhancer-driven GRNs with SCENIC+
- single-cell/trajectory-inference - Trajectory analysis for cell fate context
- single-cell/perturb-seq - Experimental perturbation data analysis with Pertpy

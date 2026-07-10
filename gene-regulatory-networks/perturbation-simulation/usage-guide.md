# Perturbation Simulation - Usage Guide

## Overview

Simulate transcription factor perturbation effects on cell state with CellOracle (and Dynamo for velocity-based fields; GEARS/CPA for response prediction). CellOracle combines a base GRN from chromatin accessibility with scRNA-seq to predict how cell identities shift under TF knockout or overexpression. The caveats that limit interpretation: these methods predict a DIRECTION, not a calibrated magnitude; they are a local linear approximation valid only near the observed manifold for a few steps; they inherit every error of the underlying GRN (CellOracle) or velocity field (Dynamo); and for perturbation-response prediction, deep models often fail to beat trivial mean/additive baselines (Ahlmann-Eltze 2025). Use them for hypothesis generation, validated against real perturbation data.

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

### Response Prediction with Baselines
> "Predict the response to an unseen double perturbation with GEARS, and compare it against the additive baseline."

## What the Agent Will Do

1. Scan accessible regions for TF binding motifs to build the base GRN
2. Import scRNA-seq data and base GRN into CellOracle
3. Fit GRN models per cell type using regularized regression
4. Simulate TF perturbation (knockout or overexpression)
5. Calculate predicted cell state shifts on the embedding
6. Visualize results as quiver plots and gradient plots

## Tips

- Direction, not magnitude - the deliverable is which way cells move; predicted amplitudes are uncalibrated. Treat outputs as hypotheses, not quantitative knockout transcriptomes.
- Local validity only - keep perturbations small and near observed states; CellOracle moves probability mass among existing states and cannot invent a new cell type. Do not crank n_propagation (default 3) for "long-range" effects.
- Inherits the substrate's errors - CellOracle is only as good as its GRN; Dynamo only as good as its RNA velocity (which is itself error-prone). Sanity-check the underlying network/field.
- Score against a development field - the CellOracle perturbation score is the inner product with the differentiation flow; it is meaningless without a defined vector field.
- Beat the baselines - for response prediction, report the mean baseline (unseen singles) and additive baseline (combinations); a model that does not beat them adds no value.
- Validate experimentally - compare predictions against Perturb-seq or published KO phenotypes; report how many predictions were tested and which failed.
- Base GRN flexibility - it can come from scATAC, bulk ATAC, or a prebuilt CellOracle base GRN; paired multiome is not required (build it in multiomics-grn).

## Related Skills

- multiomics-grn - build the CellOracle base GRN from accessibility
- scenic-regulons - regulon activity as an alternative TF-driver readout
- single-cell/trajectory-inference - pseudotime/development flow for the perturbation score
- single-cell/perturb-seq - experimental Perturb-seq, the interventional ground truth

# Batch Design Usage Guide

## Overview

This guide covers designing genomics experiments so that technical nuisance variation (batch, lane, plate, flow cell, operator, reagent lot, processing day) is balanced against the biological variable of interest and therefore estimable rather than confounded. The governing principle is that no post-hoc correction recovers a design in which batch is aliased with the biological variable, and that even on unbalanced designs aggressive correction can inject false signal. The skill covers constrained sample-to-batch assignment, the causal-graph view of what to adjust for, detection of hidden batches with surrogate variable analysis, and a decision table for downstream correction whose execution is deferred to the differential-expression and single-cell categories.

## Prerequisites

```r
# R/Bioconductor
install.packages('BiocManager')
BiocManager::install(c('sva', 'RUVSeq', 'limma', 'edgeR', 'OSAT'))

# Design optimization (CRAN)
install.packages('designit')
```

## Quick Start

Tell your AI agent what you want to do:
- "Assign 24 samples to 3 sequencing batches without confounding condition with batch"
- "Is my batch layout balanced, or is batch aliased with my condition?"
- "Can batch correction rescue my design, or do I need to redesign?"
- "Estimate how many hidden batches are in my expression data with SVA"
- "Which batch-correction method fits my design, and where do I run it?"

## Example Prompts

### Experimental Design

> "I have 12 treated and 12 control samples and 3 sequencing batches. Assign them so batch is orthogonal to condition and sex is balanced too."

> "My plate holds 16 wells but I have 4 conditions x 6 replicates. How do I split across plates without confounding?"

### Confounding and Salvage

> "All my tumor samples were sequenced in one run and all normals in another. Can ComBat fix this?"

> "My scRNA-seq has one donor per lane. How do I break the donor-versus-lane confound in the next experiment?"

### Detection and Correction Choice

> "Check my RNA-seq for hidden batch structure with SVA and tell me how to include it in the model."

> "Should I use ComBat-seq, RUVSeq, or keep batch as a covariate for my balanced design?"

## What the Agent Will Do

1. Assess the design for confounding between batch and the biological variable.
2. Produce a balanced, constrained sample-to-batch/lane/plate assignment.
3. Estimate hidden batch structure with surrogate variable analysis when warranted.
4. Recommend a correction strategy (preferring batch-in-the-model) and point to where it is executed.
5. Specify the reproducibility metadata to record so technical factors remain modelable.

## Tips

- Never run all of one condition in a single batch; balance every condition across all batches.
- Keep batch in the analysis model for inference; reserve a batch-"cleaned" matrix for visualization only.
- On unbalanced designs, aggressive correction can inject false positives; redesign rather than rely on correction.
- Adjust for confounders (common causes), not for mediators or colliders.
- Record date, reagent lot, operator, lane, and plate position for every sample.
- For scRNA-seq, pool donors per lane and demultiplex to remove the donor-versus-lane confound.

## Related Skills

- randomization-blocking - The experimental unit, randomization, and blocking concepts behind a good batch layout
- power-analysis - Account for blocking/batch factors in the power calculation
- sample-size - Balanced designs assume equal n per group
- differential-expression/batch-correction - Executes ComBat-seq/RUVSeq/SVA on real data
- single-cell/batch-integration - scRNA-seq integration (Harmony, scVI, Seurat anchors)
- machine-learning/model-validation - Batch confounding is a data-leakage source for ML

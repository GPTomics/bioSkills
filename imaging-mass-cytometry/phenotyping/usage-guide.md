# Phenotyping - Usage Guide

## Overview

Assigns cell-type identities to segmented IMC/MIBI cells from marker expression. The core caveat: a "cell type" in imaging is an inference conditioned on a segmentation guess, so the most dangerous phenotypes (CD3+CD20+ double-positives) are usually boundary artifacts, not biology -- clustering invents them as named populations while marker-dictionary classifiers refuse them, which is why the choice of method matters as much as the parameters.

## Prerequisites

```bash
pip install scanpy anndata astir scikit-learn squidpy numpy
# R FlowSOM (optional, clustering)
# BiocManager::install('FlowSOM')
```

## Quick Start

Tell your AI agent what you want to do:
- "Classify my cells with a marker dictionary and flag ambiguous cells as Unknown"
- "Cluster my cells on lineage markers and annotate the clusters"
- "Tell me whether my CD3+CD20+ population is real or a segmentation artifact"
- "Separate lineage assignment from activation-state profiling"
- "Transfer cell-type labels from my annotated reference to new samples"

## Example Prompts

### Marker-based classification
> "I can write down which markers define which cell type. Use Astir so ambiguous cells become Unknown instead of a forced call, and tell me the Unknown rate."

### Clustering
> "Cluster these cells with Leiden on lineage markers only, then show me a marker heatmap per cluster. Keep Ki67 and PD-1 out of the clustering."

### Artifact diagnosis
> "I found a CD3+CD68+ cluster. Map it back onto the image and tell me if those cells sit on T-cell/macrophage borders."

### Method choice
> "I have expert labels and bad segmentation in dense lymphoid tissue. Which phenotyping method rejects spillover double-positives best?"

## What the Agent Will Do

1. Transform the matrix with the IMC arcsinh cofactor (~1), keeping raw counts and treating zeros as real low counts.
2. Pick clustering vs classification from the decision tree (Astir for rule-based lineage; CellSighter when labels exist and spillover is heavy; Pixie/Leiden for exploration).
3. Cluster or classify on lineage markers, separating activation-state markers into a second layer.
4. Diagnose any biologically-impossible double-positive by mapping it onto the image and checking border/donor-adjacency.
5. Hand off cross-condition comparison of the resulting types to differential-analysis (patient as the unit).

## Tips

- Do not reuse the suspension-CyTOF cofactor 5; IMC single-cell means use ~1.
- A tidy double-positive cluster is a segmentation/spillover diagnosis until the image proves otherwise.
- Channel compensation (CATALYST) and lateral-spillover compensation (REDSEA) are different fixes; running one does not handle the other, and neither fixes a merged segment.
- A near-zero Astir Unknown rate is suspicious (over-confident); a high rate flags a mis-specified panel/dictionary.
- Validating clusters with the same markers used to cluster them confirms nothing -- use held-out evidence.
- The million-cell count is a red herring for any group comparison; the replicate is the patient.

## Related Skills

- cell-segmentation - double-positives diagnose the segmentation/spillover that phenotyping inherits
- data-preprocessing - arcsinh cofactor and channel spillover compensation
- spatial-analysis - phenotype labels feed neighborhood and niche analysis
- differential-analysis - comparing cell-type proportions across conditions at the patient level
- interactive-annotation - mapping clusters back onto tissue to confirm they are real
- flow-cytometry/clustering-phenotyping - FlowSOM/PhenoGraph background for suspension data

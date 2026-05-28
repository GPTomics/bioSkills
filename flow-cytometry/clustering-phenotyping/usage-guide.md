# Clustering and Phenotyping - Usage Guide

## Overview
Unsupervised clustering identifies cell populations without predefined gates - essential for high-dimensional CyTOF and spectral flow. This skill covers algorithm choice (FlowSOM is the benchmarked default; PhenoGraph and X-shift for graph/density discovery), the two rules that carry most of the correctness (cluster on type/lineage markers only, and treat UMAP/tSNE as visualization not measurement), the over-provision-then-metacluster discipline, seed dependence, and annotation by median-marker heatmap. Population definitions are made in high-dimensional space; the embedding only colors them.

## Prerequisites
```bash
# R/Bioconductor
BiocManager::install(c('CATALYST', 'FlowSOM', 'flowCore'))
# PhenoGraph (GitHub only)
# remotes::install_github('JinmiaoChenLab/Rphenograph')
```

## Quick Start
Tell your AI agent what you want to do:
- "Cluster my CyTOF data with FlowSOM and annotate the populations"
- "Run PhenoGraph and tell me how many populations it finds"
- "Make a UMAP colored by cluster and a median-marker heatmap"
- "Help me choose the number of metaclusters"

## Example Prompts
### Clustering
> "Run the CATALYST FlowSOM pipeline on my CyTOF SCE using only the lineage markers, over-provision the grid, and cap metaclusters at 20 - and set a seed so it's reproducible."
> "Cluster with PhenoGraph at k=30 and compare the number of communities to my FlowSOM metaclusters."

### Annotation and merging
> "Show me a median-marker heatmap per cluster and propose cell-type labels, then merge the over-split clusters with a curated table."
> "Which of my markers are type vs state, and why should the state markers be excluded from clustering?"

### Visualization and stability
> "Make a UMAP colored by cluster, but quantify abundances in marker space, not on the embedding."
> "This rare population only shows up at one seed - help me decide if it's real."

## What the Agent Will Do
1. Build the SCE with `prepData` and confirm panel `marker_class` (type vs state).
2. Cluster on type markers only (FlowSOM via CATALYST), over-provisioning the grid and setting a seed.
3. Annotate clusters from a median-marker heatmap and merge over-split clusters.
4. Generate a UMAP/tSNE for visualization, coloring by cluster.
5. Flag seed-unstable clusters and recommend stability/orthogonal checks before calling them novel.

## Tips
- Cluster on type/lineage markers; withhold state/functional markers for differential-state testing.
- Over-provision the SOM grid then metacluster down - clustering can merge but never split a fused node.
- Set and report the seed (`seed=` on `cluster()`); FlowSOM/PhenoGraph are stochastic.
- CATALYST `cluster()` caps metaclusters at `maxK=20` by default - raise it if you expect more.
- UMAP/tSNE distances, cluster sizes, and densities are not quantitative - never gate on the embedding.
- Use median (not mean) per-cluster expression - robust to residual doublets/spillover.
- Compensate and transform before clustering, or spillover and scale dominate distances.
- A "novel population" needs seed stability + a genuinely bimodal marker + ideally orthogonal/independent-cohort validation.

## Related Skills
- compensation-transformation - Compensate/transform before clustering
- gating-analysis - Supervised alternative; needed for rare populations
- differential-analysis - Test abundance/state of clusters between conditions
- cytometry-qc - Cluster only QC-passed events
- single-cell/clustering - Leiden/Louvain on scRNA-seq (shared graph-clustering ideas)
- imaging-mass-cytometry/phenotyping - Same CATALYST/FlowSOM conventions for imaging

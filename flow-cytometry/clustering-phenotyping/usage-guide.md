# Clustering and Phenotyping - Usage Guide

## Overview
Unsupervised clustering identifies cell populations without predefined gates. Useful for discovery and high-dimensional CyTOF/spectral flow data.

## Prerequisites
```bash
# R/Bioconductor
BiocManager::install(c('FlowSOM', 'CATALYST', 'Rphenograph'))
```

## Quick Start
Tell your AI agent what you want to do:
- "Cluster my CyTOF data with FlowSOM"
- "Identify cell populations using Phenograph"
- "Annotate clusters based on marker expression"

## Example Prompts
### Clustering
> "Run FlowSOM clustering with 20 metaclusters on my CyTOF data"
> "Apply Phenograph clustering to identify cell populations"
> "Use Leiden clustering via CATALYST"

### Visualization
> "Create a UMAP colored by FlowSOM clusters"
> "Show a heatmap of marker expression per cluster"
> "Generate a cluster abundance plot per sample"

### Annotation
> "Annotate clusters based on canonical marker expression"
> "Merge similar clusters that express the same markers"
> "Create a panel file with type and state marker annotations"

## What the Agent Will Do
1. Prepare SingleCellExperiment or flowSet for clustering
2. Select phenotyping markers (exclude state markers)
3. Apply clustering algorithm (FlowSOM, Phenograph, or Leiden)
4. Generate UMAP embedding
5. Visualize and annotate clusters

## Tips
- FlowSOM is fast and scalable; Phenograph finds rare populations
- Use type markers (lineage) for clustering, state markers for downstream analysis
- Start with more clusters than expected, then merge
- Check delta area plot to help choose number of clusters
- CATALYST provides a streamlined workflow for CyTOF

## Methods Comparison

| Method | Speed | Auto K | Best For |
|--------|-------|--------|----------|
| FlowSOM | Fast | No | Large datasets, defined biology |
| Phenograph | Slow | Yes | Discovery, rare populations |
| Leiden | Medium | Yes | Consistent community detection |

## CATALYST Panel Format
```
fcs_colname,antigen,marker_class
Yb176Di,CD45,type
Er168Di,CD3,type
Nd142Di,Ki67,state
```

## References
- FlowSOM: doi:10.1002/cyto.a.22625
- Phenograph: doi:10.1016/j.cell.2015.05.047
- CATALYST: doi:10.1101/218826

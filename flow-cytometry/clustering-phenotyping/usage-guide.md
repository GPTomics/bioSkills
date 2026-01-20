# Clustering and Phenotyping Usage Guide

## Overview

Unsupervised clustering identifies cell populations without predefined gates. Useful for discovery and high-dimensional data.

## Methods

### FlowSOM
- Self-organizing map + hierarchical clustering
- Fast, scalable
- Requires specifying # clusters

### Phenograph
- KNN graph + Louvain clustering
- Automatic # clusters
- Good for rare populations

### Leiden
- Improved Louvain
- More consistent results
- Implemented in CATALYST

## Workflow

1. **Preprocess**: Transform data (arcsinh for CyTOF)
2. **Select markers**: Use phenotyping markers
3. **Cluster**: FlowSOM or Phenograph
4. **Visualize**: UMAP colored by cluster
5. **Annotate**: Name clusters based on markers

## Choosing K (Number of Clusters)

- Start with overestimate (20-30)
- Check delta area plot
- Merge similar clusters
- Match to known biology

## Marker Selection

### Type markers
- Cell lineage definition
- Used for clustering

### State markers
- Functional states
- Used for differential analysis

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

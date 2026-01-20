# Trajectory Inference - Usage Guide

## Overview

Trajectory inference reconstructs developmental paths and orders cells by pseudotime to understand differentiation, cell state transitions, and temporal gene expression dynamics.

## Prerequisites

```r
# R packages
install.packages('BiocManager')
BiocManager::install(c('monocle3', 'slingshot', 'tradeSeq'))
```

```bash
# Python packages
pip install scvelo scanpy

# Velocyto for spliced/unspliced counts
pip install velocyto
```

## Quick Start

Tell your AI agent what you want to do:
- "Run trajectory analysis on my Seurat object"
- "Compute RNA velocity to find differentiation direction"
- "Order cells by pseudotime from stem cells"

## Example Prompts

### Trajectory Construction
> "Run Monocle3 trajectory analysis on my data"
> "Use Slingshot to find lineage curves"
> "Build a PAGA graph to find cell state transitions"

### Pseudotime
> "Compute pseudotime starting from cluster 0"
> "Show pseudotime values on the UMAP"
> "Which cells are at the beginning vs end of differentiation?"

### RNA Velocity
> "Run scVelo RNA velocity analysis"
> "Show velocity arrows on my UMAP"
> "Find driver genes for this trajectory"

### Gene Dynamics
> "Find genes that change along pseudotime"
> "Show expression of GATA1 along the trajectory"
> "Run tradeSeq to find trajectory-associated genes"

### Branch Analysis
> "Identify branch points in the trajectory"
> "Compare gene expression between the two branches"
> "Which genes drive the decision at the branch point?"

## What the Agent Will Do

1. Verify data is clustered and annotated
2. Identify starting population (root)
3. Learn trajectory graph structure
4. Compute pseudotime ordering
5. Visualize trajectory on UMAP
6. Identify genes changing along trajectory
7. Analyze branch points if present

## Tool Selection

| Tool | Best For | Language | Input |
|------|----------|----------|-------|
| Monocle3 | Complex branching trajectories | R | Seurat/SCE |
| Slingshot | Smooth lineage curves | R | SCE with clusters |
| scVelo | RNA velocity, directionality | Python | Spliced/unspliced counts |
| PAGA | Graph abstraction | Python | Scanpy AnnData |

## Tips

- **Root selection is critical** - identifies starting point for pseudotime
- **Cluster well first** - trajectory methods need good clustering
- **RNA velocity adds directionality** - use for biological direction confirmation
- **Validate with known markers** - check that marker dynamics match expectations
- **Branch points need validation** - biological knowledge helps interpret splits
- **scVelo needs spliced/unspliced** - run velocyto on BAM files first

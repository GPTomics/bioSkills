# Trajectory Inference Usage Guide

## Overview

Trajectory inference reconstructs developmental paths and orders cells by pseudotime. Choose your tool based on data type and analysis goals.

## Tool Selection

| Tool | Best For | Language | Input |
|------|----------|----------|-------|
| Monocle3 | Complex branching trajectories | R | Seurat/SCE |
| Slingshot | Smooth lineage curves | R | SCE with clusters |
| scVelo | RNA velocity, directionality | Python | Spliced/unspliced counts |
| PAGA | Graph abstraction | Python | Scanpy AnnData |

## Quick Start Prompts

- "Run trajectory analysis on my Seurat object to find differentiation paths"
- "Compute RNA velocity to determine trajectory direction"
- "Find genes that change along the pseudotime trajectory"
- "Identify branch points in my developmental trajectory"

## Workflow

1. **Cluster cells first** - Trajectory methods need good clustering
2. **Choose root** - Identify starting population (stem/progenitor)
3. **Infer trajectory** - Learn graph structure
4. **Order cells** - Compute pseudotime from root
5. **Validate** - Check marker gene dynamics

## Requirements

```bash
# R packages
install.packages('BiocManager')
BiocManager::install(c('monocle3', 'slingshot', 'tradeSeq'))

# Python packages
pip install scvelo scanpy

# Velocyto for spliced/unspliced counts
pip install velocyto
```

## Key Considerations

- **Root selection**: Critical for pseudotime interpretation
- **Branching**: Monocle3/Slingshot handle multiple lineages
- **Directionality**: RNA velocity provides biological direction
- **Validation**: Confirm with known marker dynamics

## Related Skills

- **single-cell/clustering-annotation** - Prerequisite for trajectory
- **single-cell/cell-communication** - Signaling along trajectory

# Spatial Neighbor Graphs - Usage Guide

## Overview

This skill covers building spatial neighbor graphs for spatial transcriptomics analysis using Squidpy. Neighbor graphs define which spots/cells are considered spatial neighbors for downstream analyses.

## Prerequisites

```bash
pip install squidpy scanpy
```

## Quick Start

Tell your AI agent what you want to do:
- "Build a spatial neighbor graph"
- "Create a KNN graph with 6 neighbors"

## Example Prompts

### Basic Graph
> "Build a spatial neighbor graph for my Visium data"

> "Create a Delaunay triangulation graph"

### Specific Parameters
> "Build a KNN graph with 10 neighbors"

> "Connect spots within 100 pixels of each other"

## What the Agent Will Do

1. Extract spatial coordinates
2. Build neighbor graph using specified method
3. Store connectivities and distances in adata.obsp
4. Report graph statistics

## Tips

- **Visium grid** - Use `coord_type='grid'` and `n_rings` for hexagonal grids
- **n_neighs** - 6 is typical for Visium; adjust for single-cell resolution data
- **Delaunay** - Creates natural tessellation, good for variable density
- **Radius** - Use when you want fixed distance threshold

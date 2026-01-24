# Spatial Cell-Cell Communication - Usage Guide

## Overview

This skill covers analyzing cell-cell communication in spatial transcriptomics data using ligand-receptor analysis with Squidpy. Identify which cell types are communicating and through which signaling pathways.

## Prerequisites

```bash
pip install squidpy scanpy pandas networkx matplotlib
```

## Quick Start

Tell your AI agent what you want to do:
- "Run ligand-receptor analysis on my spatial data"
- "Find which cell types are communicating"

## Example Prompts

### Basic Analysis
> "Analyze cell-cell communication in my Visium data"

> "Run ligand-receptor analysis between cell types"

### Specific Interactions
> "Find interactions between T cells and macrophages"

> "Check if CCL2-CCR2 signaling is active"

### Visualization
> "Show a heatmap of cell-cell interactions"

> "Plot the communication network"

## What the Agent Will Do

1. Verify cell type annotations exist
2. Build spatial neighbor graph
3. Run ligand-receptor permutation test
4. Filter significant interactions
5. Visualize results

## Tips

- **Cell type annotations** - Required before running communication analysis
- **Spatial neighbors** - Analysis considers spatial proximity
- **Permutations** - More permutations = more robust p-values but slower
- **Custom databases** - Can use custom ligand-receptor pair lists

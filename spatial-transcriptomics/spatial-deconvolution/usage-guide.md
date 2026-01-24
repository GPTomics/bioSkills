# Spatial Deconvolution - Usage Guide

## Overview

This skill covers estimating cell type composition in spatial transcriptomics spots using reference-based deconvolution methods like cell2location, Tangram, and RCTD.

## Prerequisites

```bash
# cell2location
pip install cell2location

# Tangram
pip install tangram-sc

# RCTD (R package)
# install.packages('spacexr')
```

## Quick Start

Tell your AI agent what you want to do:
- "Deconvolve my Visium data using this scRNA-seq reference"
- "Estimate cell type proportions in each spot"

## Example Prompts

### Basic Deconvolution
> "Run cell2location on my spatial data"

> "Use Tangram to map cell types to spatial spots"

### Visualization
> "Plot cell type proportions spatially"

> "Show the dominant cell type in each spot"

### Comparison
> "Compare deconvolution results with marker gene expression"

## What the Agent Will Do

1. Load spatial and reference scRNA-seq data
2. Find shared genes and preprocess
3. Train deconvolution model
4. Estimate cell type proportions per spot
5. Visualize results

## Tips

- **Reference quality** - Better reference = better deconvolution
- **Marker genes** - Help Tangram; cell2location learns automatically
- **N_cells_per_location** - Adjust based on tissue/platform
- **Validation** - Check correlation with known marker genes

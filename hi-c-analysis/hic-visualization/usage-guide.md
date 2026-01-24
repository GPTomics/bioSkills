# Hi-C Visualization - Usage Guide

## Overview

This skill covers visualizing Hi-C contact matrices, TADs, loops, and other genomic features using matplotlib, cooltools, and HiCExplorer.

## Prerequisites

```bash
pip install cooler cooltools matplotlib numpy
# For HiCExplorer:
conda install -c bioconda hicexplorer
```

## Quick Start

Tell your AI agent what you want to do:
- "Plot my Hi-C contact matrix"
- "Show Hi-C with TADs and loops"

## Example Prompts

### Basic Plots
> "Plot the contact matrix for chr1:50-60Mb"

> "Create a triangle plot of my Hi-C data"

### With Annotations
> "Plot Hi-C with TAD boundaries"

> "Show loops on the contact matrix"

### Comparisons
> "Compare Hi-C between treatment and control"

> "Create a split view of two samples"

## What the Agent Will Do

1. Load cooler file
2. Extract matrix for requested region
3. Create matplotlib figure
4. Add annotations (TADs, loops) if requested
5. Save figure

## Tips

- **Log scale** - Use LogNorm for contact matrices
- **Color limits** - vmin/vmax control dynamic range
- **Resolution** - Higher resolution = more detail but slower
- **Triangle plots** - Better for linear arrangement with tracks

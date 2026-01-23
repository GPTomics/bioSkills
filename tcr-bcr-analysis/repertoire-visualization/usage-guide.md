# Repertoire Visualization - Usage Guide

## Overview

Create publication-quality visualizations of immune repertoire data including V-J circos plots, clone tracking, diversity comparisons, and clonotype networks.

## Prerequisites

```bash
pip install matplotlib seaborn networkx python-Levenshtein
# or
conda install -c conda-forge matplotlib seaborn networkx python-levenshtein
```

## Quick Start

Tell your AI agent:
- "Create a circos plot of V-J gene usage"
- "Track my top clones across timepoints"
- "Plot diversity comparison between conditions"
- "Generate a repertoire overlap heatmap"

## Example Prompts

### Gene Usage

> "Create a V-J pairing circos plot"

> "Plot V gene usage by sample"

> "Compare J gene frequencies between groups"

### Clone Dynamics

> "Track the top 10 clones over time"

> "Show clonal expansion in my samples"

> "Plot frequency change of shared clonotypes"

### Diversity

> "Compare Shannon diversity between treatment and control"

> "Create a boxplot of Gini coefficients by condition"

> "Plot rarefaction curves for my samples"

### Comparison

> "Generate overlap heatmap for all sample pairs"

> "Create a network of similar clonotypes"

> "Show public vs private clonotypes"

## What the Agent Will Do

1. Load processed clonotype data
2. Calculate metrics for visualization
3. Create appropriate plot type
4. Apply formatting for publication quality
5. Save in requested format (PDF, PNG)

## Tips

- **Circos plots** work best with aggregated V/J gene families
- **Clone tracking** - limit to top 10-20 clones for readability
- **Heatmaps** - hierarchical clustering helps identify patterns
- **Networks** - filter by similarity threshold to reduce complexity
- **Color** - use colorblind-friendly palettes

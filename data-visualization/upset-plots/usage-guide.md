# UpSet Plots - Usage Guide

## Overview

UpSet plots visualize set intersections more effectively than Venn diagrams, especially for 4+ sets. They show intersection sizes as bar charts with a matrix indicating which sets participate.

## Prerequisites

```r
# R
install.packages('UpSetR')
```

```bash
# Python
pip install upsetplot matplotlib pandas
```

## Quick Start

Tell your AI agent:
- "Create an UpSet plot of my gene set overlaps"
- "Compare peak overlaps across my ChIP-seq samples"
- "Show which genes are shared between my DE results"
- "Visualize the intersection of my variant filter sets"

## Example Prompts

### Basic UpSet Plots

> "Create an UpSet plot showing overlaps between my 5 gene lists"

> "Make an UpSet plot of peaks shared across ChIP-seq replicates"

### Customization

> "Highlight the intersection of SetA and SetB in red"

> "Sort intersections by degree instead of frequency"

> "Add set size bars to my UpSet plot"

### With Metadata

> "Add boxplots showing log fold change for each intersection"

> "Color intersections by whether they contain significant genes"

## What the Agent Will Do

1. Convert input data to appropriate format (list or binary matrix)
2. Determine number of sets and optimal layout
3. Configure sorting and display options
4. Add any requested highlights or annotations
5. Export plot in requested format

## Tips

- **Use UpSet over Venn** for 4+ sets - Venn diagrams become unreadable
- **Sort by frequency** (default) to see largest intersections first
- **Sort by degree** to group by number of participating sets
- **Queries** in UpSetR let you highlight specific intersections
- **Limit intersections** with nintersects if plot is too crowded
- **Add attributes** to show metadata distributions per intersection

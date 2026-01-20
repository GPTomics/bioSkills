# Cell-Cell Communication - Usage Guide

## Overview

Cell-cell communication analysis identifies ligand-receptor interactions between cell types from scRNA-seq data to understand tissue organization, signaling networks, and intercellular coordination.

## Prerequisites

```r
# CellChat
devtools::install_github('sqjin/CellChat')

# NicheNet
devtools::install_github('saeyslab/nichenetr')
# Download NicheNet databases from Zenodo
```

```bash
# LIANA (Python)
pip install liana
```

## Quick Start

Tell your AI agent what you want to do:
- "Find ligand-receptor interactions between cell types"
- "Compare cell communication between conditions"
- "Identify signaling pathways between macrophages and T cells"

## Example Prompts

### Interaction Discovery
> "Run CellChat to find all ligand-receptor interactions"
> "Identify which cell types communicate the most"
> "Find interactions specific to fibroblast-epithelial crosstalk"

### Pathway Analysis
> "Which signaling pathways are active between these cell types?"
> "Show a network of WNT signaling interactions"
> "Find all TGF-beta pathway communications"

### Condition Comparison
> "Compare cell communication between control and disease"
> "Which interactions are gained or lost in treatment?"
> "Show differential signaling between conditions"

### Target Prediction
> "Use NicheNet to find ligands affecting gene expression in T cells"
> "Which ligands from macrophages activate inflammatory genes?"
> "Predict downstream targets of fibroblast-derived signals"

### Visualization
> "Create a chord diagram of cell communications"
> "Show a heatmap of interaction strengths"
> "Plot the ligand-receptor network"

## What the Agent Will Do

1. Load annotated scRNA-seq data with cell types
2. Select appropriate ligand-receptor database
3. Calculate communication probabilities
4. Identify significant interactions
5. Aggregate at pathway level
6. Visualize networks and interactions
7. Compare conditions if applicable

## Tool Selection

| Tool | Strengths | Language | Best For |
|------|-----------|----------|----------|
| CellChat | Comprehensive, pathway-level | R | General CCC, comparison |
| NicheNet | Ligand-target prediction | R | Functional impact on receivers |
| LIANA | Multiple methods, consensus | Python | Robust rankings, multi-sample |

## Tips

- **Cell type annotation quality** directly affects results
- **Expression thresholds** filter noise but may miss real interactions
- **Database choice** affects coverage (CellChatDB vs OmniPath vs custom)
- **Multiple testing** - adjust p-values for many comparisons
- **Validate top hits** - check ligand and receptor expression patterns
- **Spatial context helps** - pair with spatial data if available

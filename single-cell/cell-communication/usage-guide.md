# Cell-Cell Communication Usage Guide

## Overview

Cell-cell communication analysis identifies ligand-receptor interactions between cell types from scRNA-seq data to understand tissue organization and signaling networks.

## Tool Selection

| Tool | Strengths | Language | Best For |
|------|-----------|----------|----------|
| CellChat | Comprehensive, pathway-level | R | General CCC, comparison |
| NicheNet | Ligand-target prediction | R | Functional impact on receivers |
| LIANA | Multiple methods, consensus | Python | Robust rankings, multi-sample |

## Quick Start Prompts

- "Identify ligand-receptor interactions between macrophages and T cells"
- "Compare cell communication between control and treatment conditions"
- "Find which ligands from fibroblasts activate genes in epithelial cells"
- "Run LIANA to get consensus cell-cell communication scores"

## Workflow

1. **Prepare data** - Annotated scRNA-seq with cell types
2. **Select database** - CellChatDB, NicheNet resources, or consensus
3. **Identify interactions** - Compute communication probabilities
4. **Visualize** - Network plots, chord diagrams, heatmaps
5. **Compare** - Differential communication between conditions

## Requirements

```r
# CellChat
devtools::install_github('sqjin/CellChat')

# NicheNet
devtools::install_github('saeyslab/nichenetr')
# Download NicheNet databases from Zenodo
```

```bash
# LIANA
pip install liana
```

## Key Considerations

- **Cell type annotation quality** directly affects results
- **Expression thresholds** filter noise vs lose real interactions
- **Database choice** affects coverage and specificity
- **Multiple testing** - many interactions, adjust p-values

## Related Skills

- **single-cell/clustering-annotation** - Prerequisite cell types
- **spatial-transcriptomics/spatial-communication** - Adds spatial context

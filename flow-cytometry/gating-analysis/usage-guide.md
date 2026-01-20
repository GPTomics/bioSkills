# Gating Analysis - Usage Guide

## Overview
Gating defines cell populations by drawing boundaries in parameter space. Essential for identifying cell subsets in flow cytometry and CyTOF data.

## Prerequisites
```bash
# R/Bioconductor
BiocManager::install(c('flowCore', 'flowWorkspace', 'openCyto', 'flowDensity'))
```

## Quick Start
Tell your AI agent what you want to do:
- "Gate lymphocytes from my flow data"
- "Apply an automated gating strategy"
- "Extract CD4+ T cells using a standard gating hierarchy"

## Example Prompts
### Manual Gating
> "Create a polygon gate for lymphocytes based on FSC-A vs SSC-A"
> "Apply a rectangular gate for CD3+ cells"
> "Set up a quadrant gate for CD4 vs CD8"

### Automated Gating
> "Use flowDensity to automatically gate the major populations"
> "Apply an openCyto gating template to my flowSet"
> "Find the optimal gate boundaries using density-based clustering"

### Gating Hierarchies
> "Apply a standard T cell gating strategy: singlets -> live -> CD3+ -> CD4/CD8"
> "Extract population statistics from my gating hierarchy"
> "Export gated populations as new FCS files"

## What the Agent Will Do
1. Load gating packages (flowWorkspace, openCyto, or flowDensity)
2. Create GatingSet from flowSet
3. Apply gates (manual or automated) in hierarchical order
4. Visualize gating strategy
5. Extract population statistics or gated cells

## Tips
- Always gate singlets first (FSC-A vs FSC-H)
- Rectangular gates are fast and reproducible
- Polygon gates offer more flexibility for complex populations
- openCyto templates enable reproducible automated gating
- flowDensity uses density peaks for data-driven boundaries

## Gate Types

| Type | Use Case |
|------|----------|
| Rectangular | Simple threshold on 1-2 parameters |
| Polygon | Arbitrary shapes for complex populations |
| Quadrant | 4 populations from 2 parameters (e.g., CD4/CD8) |
| Boolean | Combine gates with AND, OR, NOT |

## Automated Gating Tools

| Tool | Approach |
|------|----------|
| flowDensity | Density-based automatic gating |
| openCyto | Template-based reproducible gating |
| flowClust | Model-based (mixture models) |

## References
- flowWorkspace: doi:10.1186/s12859-018-2425-9
- openCyto: doi:10.1371/journal.pcbi.1003806

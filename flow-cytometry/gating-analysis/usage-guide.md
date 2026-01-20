# Gating Analysis Usage Guide

## Overview

Gating defines cell populations by drawing boundaries in parameter space. Essential for identifying cell subsets.

## Gate Types

### Rectangular
- Simple threshold on 1-2 parameters
- Fast, reproducible
- May miss oblique populations

### Polygon
- Arbitrary shape
- More flexible
- Harder to reproduce

### Quadrant
- 4 populations from 2 parameters
- Common for CD4/CD8

### Boolean
- Combine gates with AND, OR, NOT
- Complex population definitions

## Gating Strategies

### Standard T Cell Panel
1. FSC-A vs SSC-A → Lymphocytes
2. FSC-A vs FSC-H → Singlets
3. CD45 → Leukocytes
4. CD3 → T cells
5. CD4 vs CD8 → Subsets

## Automated vs Manual

### Manual
- Expert-defined
- Consistent criteria
- Time-consuming

### Automated
- Data-driven
- Reproducible
- May need validation

## Tools

| Tool | Approach |
|------|----------|
| flowDensity | Density-based |
| openCyto | Template-based |
| flowClust | Model-based |

## References

- flowWorkspace: doi:10.1186/s12859-018-2425-9
- openCyto: doi:10.1371/journal.pcbi.1003806

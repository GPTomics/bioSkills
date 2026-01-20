# Phenotyping Usage Guide

## Overview

Phenotyping assigns cell type identities based on marker expression. Critical for biological interpretation of IMC data.

## Approaches

### Manual Gating
- Define thresholds for each marker
- Boolean logic for cell types
- Most interpretable, requires expertise

### Clustering
- Unsupervised grouping
- PCA + Leiden/Louvain
- FlowSOM for cytometry-style

### Automated
- Reference-based classification
- Transfer learning from labeled data
- Faster, requires good reference

## Marker Panel Design

### Lineage Markers
- Pan-immune: CD45
- T cells: CD3, CD4, CD8
- B cells: CD20
- Macrophages: CD68, CD163
- Epithelial: E-cadherin, pan-CK

### Functional Markers
- Activation: Ki67, PD-1, PD-L1
- Cytotoxicity: Granzyme B

## Best Practices

1. **Transform data** before clustering (arcsinh)
2. **Select relevant markers** for phenotyping
3. **Validate** with known biology
4. **Report confidence** of assignments

## Common Cell Types

| Type | Key Markers |
|------|-------------|
| CD8 T | CD45+CD3+CD8+ |
| CD4 T | CD45+CD3+CD4+ |
| B cell | CD45+CD20+ |
| Macrophage | CD45+CD68+ |
| Epithelial | E-cadherin+ |

## References

- FlowSOM: doi:10.1002/cyto.a.22625
- Phenograph: doi:10.1016/j.cell.2015.05.047

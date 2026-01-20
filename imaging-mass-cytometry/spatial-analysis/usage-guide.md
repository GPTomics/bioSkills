# Spatial Analysis Usage Guide

## Overview

Spatial analysis reveals tissue organization by analyzing where cell types are located relative to each other.

## Key Analyses

### Neighborhood Enrichment
- Are cell type pairs found together more than expected?
- Permutation-based z-scores
- Positive: attraction, Negative: avoidance

### Co-occurrence
- How does co-occurrence change with distance?
- Reveals interaction ranges

### Ripley's Statistics
- Tests for clustering vs random distribution
- L(r) function commonly used

## Graph Construction

### Delaunay Triangulation
- Connects neighboring cells
- No distance threshold needed
- Good for densely packed tissues

### Radius-Based
- Connect cells within distance r
- Choose r based on cell size (~2-3 cell diameters)
- Good for sparse tissues

## Interpretation

| Z-score | Meaning |
|---------|---------|
| > 2 | Significant enrichment (attraction) |
| < -2 | Significant depletion (avoidance) |
| -2 to 2 | Random association |

## Best Practices

1. **Match radius** to biological question
2. **Consider density** (normalize for cell abundance)
3. **Test multiple radii** for robustness
4. **Report effect sizes** not just p-values

## References

- squidpy: doi:10.1038/s41592-021-01358-2
- Ripley's K: doi:10.1111/j.2517-6161.1977.tb01615.x

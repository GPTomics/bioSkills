# Compensation and Transformation Usage Guide

## Overview

Compensation removes spectral spillover between fluorochromes. Transformation normalizes dynamic range for visualization and analysis.

## Compensation

### Why?
Fluorochromes emit into multiple detectors. Spillover must be mathematically removed.

### Matrix
Square matrix with detectors on both axes. Diagonal = 1, off-diagonal = spillover coefficients.

### Best Practices
- Use single-stained controls
- Match voltages to experiment
- Verify with FMO controls

## Transformations

### Biexponential (Logicle)
- Standard for conventional flow
- Handles negative values (compensation artifacts)
- ~5 decades of display

### Arcsinh
- Standard for mass cytometry (CyTOF)
- cofactor = 5 is typical
- Formula: asinh(x/cofactor)

### Log
- Classic transformation
- Cannot handle negative/zero values
- Less common now

## When to Transform

| Data Type | Transform |
|-----------|-----------|
| Conventional flow | Logicle (biexponential) |
| CyTOF | Arcsinh (cofactor 5) |
| Spectral flow | Logicle after unmixing |

## Verification

1. Check bivariate plots post-compensation
2. Verify no spillover-driven false positives
3. Check spread (compensation increases CV)

## References

- Logicle: doi:10.1002/cyto.a.20258
- CyTOF normalization: doi:10.1002/cyto.a.22271

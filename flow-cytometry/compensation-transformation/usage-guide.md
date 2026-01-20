# Compensation and Transformation - Usage Guide

## Overview
Compensation removes spectral spillover between fluorochromes. Transformation normalizes dynamic range for visualization and analysis.

## Prerequisites
```bash
# R/Bioconductor
BiocManager::install(c('flowCore', 'flowWorkspace'))
```

## Quick Start
Tell your AI agent what you want to do:
- "Apply compensation to my flow cytometry data"
- "Transform my CyTOF data with arcsinh"
- "Calculate compensation matrix from single-stained controls"

## Example Prompts
### Compensation
> "Apply the compensation matrix from the FCS keywords to my data"
> "Build a spillover matrix from my single-stained controls"
> "Verify compensation with bivariate plots"

### Transformation
> "Apply logicle transformation to my compensated flow data"
> "Transform my CyTOF data with arcsinh cofactor 5"
> "Estimate optimal logicle parameters for each channel"

### Quality Control
> "Check my compensation with FMO controls"
> "Show before/after plots for compensation verification"

## What the Agent Will Do
1. Extract or compute compensation/spillover matrix
2. Apply compensation to correct for spectral overlap
3. Select appropriate transformation (logicle for flow, arcsinh for CyTOF)
4. Apply transformation to normalize dynamic range
5. Generate verification plots

## Tips
- Conventional flow: Use logicle (biexponential) transformation
- CyTOF: Use arcsinh with cofactor = 5
- Always verify compensation with bivariate plots
- Use single-stained controls with voltages matching the experiment
- FMO (Fluorescence Minus One) controls help verify gating boundaries

## Transformation Guide

| Data Type | Transform | Notes |
|-----------|-----------|-------|
| Conventional flow | Logicle | Handles negative values from compensation |
| CyTOF | Arcsinh (cofactor 5) | Standard for mass cytometry |
| Spectral flow | Logicle | Apply after unmixing |

## References
- Logicle: doi:10.1002/cyto.a.20258
- CyTOF normalization: doi:10.1002/cyto.a.22271

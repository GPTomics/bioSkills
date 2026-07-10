# Compensation and Transformation - Usage Guide

## Overview
Compensation removes fluorophore spillover on conventional cytometers (a linear matrix subtraction); spectral cytometers instead solve an overdetermined unmixing of the full emission spectrum, with autofluorescence modeled as an extra signature. Transformation (logicle/biexponential for fluorescence, arcsinh for mass) then stabilizes variance so populations separate. The key expert point this skill encodes: neither compensation nor unmixing removes spreading error - the Poisson-driven widening of negative populations - which is fixed at panel design, not downstream. It also covers the compensate-then-transform ordering, which cannot be swapped, and the arcsinh cofactor decision.

## Prerequisites
```bash
# R/Bioconductor
BiocManager::install(c('flowCore', 'flowStats', 'flowWorkspace', 'CATALYST'))
# Optional per-channel cofactor optimization
BiocManager::install('flowVS')
```

## Quick Start
Tell your AI agent what you want to do:
- "Apply the recorded spillover matrix and logicle-transform my fluorescence channels"
- "Estimate a compensation matrix from my single-stain controls"
- "Arcsinh-transform my CyTOF data with the right cofactor"
- "My data is from a spectral cytometer - unmix it rather than compensate"

## Example Prompts
### Compensation
> "These FACS files have a $SPILLOVER keyword from the cytometer. Apply it, then estimate and apply a logicle transform on the fluorescent channels, and show me a before/after CD3 vs CD4 plot."
> "Build a spillover matrix from my single-stained controls and unstained tube - my controls are beads, so flag if any are dimmer than the experimental sample."

### Transformation and cofactor
> "Arcsinh-transform this CyTOF data with cofactor 5, then check whether my dimmest marker is being compressed - if so, suggest a per-channel cofactor."
> "This is high-parameter fluorescence going into FlowSOM - what cofactor or logicle settings should I use, and why not cofactor 5?"

### Spectral and panel design
> "Explain whether I should compensate or unmix this Aurora spectral data, and what the unstained control is for."
> "My dim marker sits in a channel that catches spillover from a bright lineage marker - is this fixable by compensation, or a panel-design problem?"

## What the Agent Will Do
1. Determine instrument type (conventional -> compensate; spectral -> unmix; mass -> arcsinh).
2. Retrieve or estimate the spillover matrix (recorded keyword, single-stain controls, or AutoSpill).
3. Apply compensation on raw/linear data (never after a transform).
4. Choose and apply the transform (logicle via estimateLogicle for fluorescence; arcsinh cofactor 5 for CyTOF, ~150 or per-channel for fluorescence).
5. Verify with before/after bivariate plots and flag spreading-error or over-compensation artifacts.

## Tips
- Compensation corrects the mean; it cannot remove spreading error - do not hand-tune the matrix to flatten a smeared negative (that is falsification).
- Compensate BEFORE transforming; run `estimateLogicle` on compensated data so the linearization width matches the real negative spread.
- Negative values after compensation are expected and meaningful - logicle/arcsinh handle them; log cannot.
- Cofactor 5 is for CyTOF; ~150 is the fluorescence convention (not a derived optimum); optimize per-channel with flowVS when a dim marker drives a conclusion.
- Single-stain compensation controls must be at least as bright as the experimental sample and use the same fluorophore/voltages.
- Spectral unmixing requires single-stain reference spectra plus an unstained autofluorescence control; it is a different model from compensation, not "compensation with more detectors."

## Related Skills
- fcs-handling - Load FCS and retrieve the spillover keyword first
- gating-analysis - Gate on the transformed, compensated scale
- clustering-phenotyping - Cluster compensated, transformed data
- cytometry-qc - QC before and after preprocessing
- imaging-mass-cytometry/data-preprocessing - Shared arcsinh/metal-channel conventions
- spatial-transcriptomics/spatial-proteomics - Spectral/metal unmixing context

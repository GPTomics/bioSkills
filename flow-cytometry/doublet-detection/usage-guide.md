# Doublet Detection - Usage Guide

## Overview
Doublets are events where two or more cells traverse the laser (or ionize) together. They appear as false intermediate or double-positive populations and must be removed before clustering or quantification. The discriminating signal in flow is the non-proportionality of pulse Area vs Height (singlets fall on the FSC-A vs FSC-H diagonal; doublets deflect off it) plus Width - not a one-dimensional area gate. On CyTOF there is no scatter, so doublets are caught by DNA intercalator content and ion-cloud Gaussian/Event_length parameters. A key expert caveat: scatter gating is necessary but not sufficient - heterotypic conjugates survive it and masquerade as double-positive populations.

## Prerequisites
```bash
# R/Bioconductor
BiocManager::install(c('flowCore', 'flowDensity', 'CATALYST'))
```

## Quick Start
Tell your AI agent what you want to do:
- "Gate singlets using FSC-A vs FSC-H"
- "Remove doublets from all samples and report the rate per sample"
- "Identify doublets in my CyTOF data using DNA and Event_length"
- "This double-positive population looks suspicious - is it a doublet artifact?"

## Example Prompts
### Flow / spectral
> "Create a polygon singlet gate along the FSC-A vs FSC-H diagonal (not a rectangle, so I don't keep off-diagonal doublets) and apply it across my flowSet, reporting the doublet rate per sample."
> "My instrument records FSC-W - add a width-based singlet criterion on top of the A-vs-H gate."

### CyTOF
> "Gate intercalator-positive single cells using the Ir191/Ir193 DNA channels and Event_length, and explain why I can't use FSC-A here."
> "Use the Gaussian discrimination parameters to remove ion-cloud fusions that DNA content alone misses."

### Diagnosing artifacts
> "I have a CD3+CD14+ cluster sitting between my T-cell and monocyte clusters - check whether it's a heterotypic doublet by looking at parental marker intensity."

## What the Agent Will Do
1. Identify the right channels (FSC-A/H[/W] for flow; DNA + Event_length/Gaussian for CyTOF).
2. Build a diagonal singlet gate (polygon) or DNA/Gaussian filter.
3. Apply across samples and report per-sample doublet rates.
4. Flag suspicious double-positive populations as possible residual heterotypic doublets.
5. Return cleaned data for downstream clustering.

## Tips
- Use the FSC-A vs FSC-H diagonal, not a 1D area gate - doublets overlap singlets in area but not in the A-H relationship.
- Prefer a diagonal polygon to a rectangle; a box keeps off-diagonal doublets and clips large singlets.
- CyTOF has no scatter: use DNA intercalator (Ir191/193) and Gaussian/Event_length parameters.
- Scatter gating is necessary but not sufficient - heterotypic conjugates (T:monocyte) survive and look like double-positives; their lineage markers look comparable to true singlets, so the tell is an elevated shared marker (e.g. CD45) and the definitive resolver is imaging flow, not a lineage-intensity check.
- Expect ~1-5% doublets in PBMCs, higher in tissue digests; rates far above are a prep-quality flag, not a removal threshold.
- Cytometry doublet removal is gating-based; DoubletFinder/Scrublet are droplet-scRNA-seq methods with limited transfer here.

## Related Skills
- cytometry-qc - Run first: flow-rate/signal/margin cleaning
- bead-normalization - CyTOF drift correction after doublet removal
- fcs-handling - Load FCS files
- gating-analysis - Where singlet discrimination sits in the gating hierarchy
- clustering-phenotyping - Downstream analysis after doublet removal
- single-cell/doublet-detection - Droplet scRNA-seq doublet methods (different principle)

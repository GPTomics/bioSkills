# Cell-Type Deconvolution

## Overview

A bulk blood or tissue methylome is a cell-type-composition-weighted average of its constituent cell types. Because most CpG variation is between cell types rather than between conditions, a difference in composition between cases and controls reads as a difference in methylation - the single biggest confounder in epigenome-wide association studies. This skill estimates per-sample cell-type fractions from a clean beta matrix and uses them to defuse that confounder, either as EWAS covariates or by modeling which cell type drives a signal. The recurring lesson: a reference-based fraction is a projection of a sample onto cell types purified elsewhere, on another platform, in another tissue - it is only as accurate as the reference matches the sample.

## Prerequisites

Install the R/Bioconductor packages:

```r
BiocManager::install(c('EpiDISH', 'minfi', 'FlowSorted.Blood.EPIC',
                       'FlowSorted.CordBloodCombined.450k', 'TOAST'))
# 12-cell extended blood: BiocManager::install('FlowSorted.BloodExtended.EPIC')
# reference-free ReFACTor: install.packages('TCA')
# EpiSCORE (solid tissue): remotes::install_github('aet21/EpiSCORE')
```

Conceptual and data prerequisites:
- A clean, normalized beta matrix (CpGs x samples) or an `RGChannelSet` for `estimateCellCounts2`. Read-in, normalization, detection-p QC, and EPICv2 replicate-probe collapse belong to array-preprocessing.
- The reference package MUST match tissue, age, and array platform. EPIC vs 450K use different IDOL libraries; cord blood needs a reference carrying nRBC; solid tissue cannot use a flat blood reference.
- Reference experiment-data packages are large (hundreds of MB) and download on first use.
- `estimateCellCounts2` returns a `Neu` column where the older `minfi::estimateCellCounts` returns `Gran` - the label changes downstream column names.

## Quick Start

Tell your AI agent what you want to do:
- "Estimate blood cell-type fractions from my EPIC beta matrix with EpiDISH RPC"
- "Run estimateCellCounts2 with the IDOL library on my RGChannelSet"
- "Adjust my whole-blood EWAS for estimated cell composition"
- "Deconvolve my solid-tissue methylation with hierarchical EpiDISH"
- "Find which cell type drives a differential methylation signal with CellDMC"

## Example Prompts

### Reference-based blood deconvolution
> "I have a normalized EPIC beta matrix for 80 whole-blood samples. Estimate the six major immune cell-type fractions with the IDOL-optimized library and give me a samples-by-cell-type table I can use as EWAS covariates."

### Cord blood
> "These are newborn cord-blood 450K samples. Deconvolve them with a cord-blood reference that includes nucleated red blood cells, not an adult blood reference."

### Solid tissue
> "Deconvolve my tumor methylation into epithelial, fibroblast, and immune fractions, then break the immune fraction into subtypes using hierarchical EpiDISH."

### Composition as a covariate
> "Help me add the estimated cell fractions to my EWAS design matrix without introducing collinearity, and explain whether composition is a confounder or a mediator here."

### Cell-type-resolved EWAS
> "Run CellDMC to find which cell type carries my smoking-associated methylation signal, and tell me which calls I should distrust because the cell type is rare."

### Reference-free fallback
> "I have no matched reference for this tissue. Use ReFACTor to derive cell-composition covariates and check whether my top hits survive the correction."

## What the Agent Will Do

1. Confirm the input is a clean beta matrix or RGChannelSet and that the array platform (450K / EPIC / EPICv2) is known; route preprocessing and EPICv2 replicate collapse to array-preprocessing if needed.
2. Select a reference matched to tissue, age, and platform - blood IDOL (6 or 12 cell), cord-blood-with-nRBC, or a solid-tissue reference (hepidish / EpiSCORE).
3. Run reference-based deconvolution (EpiDISH RPC or estimateCellCounts2) and return a per-sample fraction table, or reference-free correction (ReFACTor / RefFreeEWAS) when no reference matches.
4. Advise on using the fractions: drop one cell type to avoid collinearity when entering them as covariates, or model a phenotype x fraction interaction (CellDMC / TCA / TOAST) to attribute a signal to a cell type.
5. Flag the failure modes: silent redistribution of missing cell types, platform-mismatched libraries, compositional collinearity, and underpowered rare-cell csDM calls.
6. Hand off composition covariates to differential-cpg-testing / ewas-design, and the IEAA cell-count adjustment to epigenetic-clocks.

## Tips

- Match the reference to tissue, age, and platform before anything else - a mismatched reference returns confident, wrong fractions with no error.
- Cord blood always needs nRBC; an adult blood reference on newborns is a classic silent error.
- Prefer RPC over the textbook-default CP in EpiDISH; it is more robust to noisy CpGs.
- When entering fractions as covariates, drop one reference cell type (or use a compositional transform) - all K fractions sum to ~1 and are collinear.
- Reference-free correction can absorb real biology; use it only when no reference exists and confirm top hits survive.
- A cell-type-specific call from bulk is a hypothesis about what to sort next, not a finding; validate the attributions your conclusions rest on in sorted or single-cell data.
- Report each cell type's mean fraction alongside any cell-type-resolved result; distrust specific calls for cell types below ~5%.

## Related Skills

- array-preprocessing - Provides the clean beta matrix deconvolution consumes
- ewas-design - Cell-fraction covariate strategy (confounder vs mediator)
- epigenetic-clocks - IEAA: adjust the clock for estimated cell composition
- differential-cpg-testing - Uses cell fractions as design-matrix covariates
- single-cell/preprocessing - scRNA atlases for reference building (EpiSCORE) and ground truth
- machine-learning/biomarker-discovery - Predictive-model boundary
- workflows/methylation-pipeline - End-to-end pipeline

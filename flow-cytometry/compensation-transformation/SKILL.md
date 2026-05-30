---
name: bio-flow-cytometry-compensation-transformation
description: Corrects fluorophore spillover (conventional compensation) or spectral overlap (spectral unmixing) and applies variance-stabilizing transforms (logicle/biexponential, arcsinh, log) for flow and mass cytometry. Covers spillover-matrix estimation from single-stain controls, AutoSpill, the spillover spreading matrix and why panel design (not compensation) bounds resolution, compensate-then-transform ordering, and arcsinh cofactor choice (5 for CyTOF, ~150 for fluorescence, per-channel via flowVS). Use when correcting spectral overlap, preparing data for gating/clustering, choosing logicle vs arcsinh, deciding a cofactor, or distinguishing compensation from spectral unmixing.
tool_type: r
primary_tool: flowCore
---

## Version Compatibility

Reference examples tested with: flowCore 2.14+, flowStats 4.14+, flowWorkspace 4.14+, CATALYST 1.26+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

Notes that bite: `estimateLogicle()` lives in flowWorkspace (not flowCore). `flowCore::spillover()` on a flowFrame returns a LIST of keyword matrices (index `[[1]]`); `flowStats::spillover()` on single-stain controls returns the matrix DIRECTLY (not a list) - do not index it with `$`.

If code throws an error, introspect the installed package and adapt rather than retrying.

# Compensation and Transformation

**"Compensate and transform my cytometry data"** -> Remove spillover (matrix subtraction, conventional) or unmix the full spectrum (least squares, spectral), then apply a transform so populations separate.
- R (conventional): `flowCore::compensate()` then `flowWorkspace::estimateLogicle()` + `flowCore::transform()`
- R (CyTOF/mass): `CATALYST::prepData(..., transform=TRUE, cofactor=5)` (arcsinh)
- R (spectral): linear UNMIXING, not compensation - see the taxonomy

## The Single Most Important Modern Insight -- Compensation Corrects the Mean; It Cannot Remove Spreading Error

Conventional compensation inverts a square spillover matrix (peak-channel subtraction); spectral cytometry solves an OVERDETERMINED least-squares unmix over all detectors, with autofluorescence modeled as an extra "fluorophore." Both correct the population MEAN. Neither removes **spreading error** - the widening of a negative population in a spillover detector that arises from the Poisson counting statistics of the spilled-in photons (Roederer 2001 *Cytometry* 45:194; Nguyen 2013 *Cytometry A* 83:306). Compensation does not INTRODUCE spreading; it makes the pre-existing variance visible by re-centering means. The corollaries are load-bearing: (1) a smeared negative cannot be fixed by tuning the matrix - over-compensating to flatten it is data falsification; (2) spreading is fixed at PANEL DESIGN (the Spillover Spreading Matrix identifies which detector pairs to avoid for co-expressed/dim markers), never downstream; (3) calling spectral unmixing "compensation" is a category error - it is a different, overdetermined model.

## Method Taxonomy

| Method | What it does | When to use | Fails when |
|--------|--------------|-------------|------------|
| Acquisition-recorded `$SPILLOVER` | applies the cytometer-computed matrix | trustworthy single-stain setup at acquisition | controls were wrong/missing |
| Computed compensation (`flowStats::spillover`) | estimates spillover from single-stain controls (medians) | conventional flow, controls available | poor/dim/contaminated controls |
| AutoSpill (Roca 2021 *Nat Commun* 12:2890) | robust-regression matrix + iterative refinement; AF as endogenous dye | high-parameter panels; messy controls | reference implementation/setup unavailable |
| Spectral unmixing (OLS/WLS/Poisson) | least-squares unmix full spectrum vs reference spectra + AF | spectral cytometers (Aurora, ID7000) | wrong/heterogeneous AF; collinear spectra |
| Logicle / biexponential | display + analysis transform, handles negatives | fluorescence flow | wrong `w` clips the negative population |
| arcsinh | variance-stabilizing transform | CyTOF/mass; computational pipelines | wrong cofactor compresses dim markers |
| log10 | legacy | rarely; strictly positive data | any negative values after compensation |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Conventional flow, `$SPILLOVER` present | apply recorded matrix -> estimateLogicle | trust acquisition controls; logicle handles negatives |
| Conventional flow, no matrix | compute via `flowStats::spillover` from single-stains (or AutoSpill) | controls drive the matrix; AutoSpill for >12 colors |
| Spectral cytometer | UNMIX (do NOT compensate), then arcsinh at ~150/per-channel (NOT 5) | overdetermined system; spectral data is fluorescence-scale, not ion counts |
| CyTOF / mass | arcsinh cofactor 5; spillover via CATALYST `compCytof` if needed | metals barely spill (~1-4%), but oxide/impurity is real |
| Dim marker driving a borderline call | test per-channel cofactor (flowVS) | a fixed cofactor can manufacture/erase the population |

## Compensate-Then-Transform Ordering (load-bearing)

Compensation/unmixing is LINEAR and must run on untransformed data; applying it after a nonlinear transform is mathematically invalid. `estimateLogicle()` must run on ALREADY-COMPENSATED data so the `w`/`a` parameters reflect the post-compensation negative spread. Negative values after compensation are expected and meaningful - do NOT clip to zero before transforming (handling negatives is the entire reason logicle/arcsinh exist; log cannot).

## Apply or Compute Compensation

**Goal:** Apply the recorded matrix, or estimate one from single-stain controls.

**Approach:** `compensate()` takes a `compensation` object built from the matrix; `flowStats::spillover()` estimates from single-stain controls and returns the matrix directly.

```r
library(flowCore)

comp <- compensation(spillover(fcs)[[1]])      # flowCore: flowFrame -> list of keyword matrices
fcs_comp <- compensate(fcs, comp)

library(flowStats)
ctrls <- read.flowSet(list.files('controls', pattern = '\\.fcs$', full.names = TRUE))
comp_matrix <- spillover(ctrls, unstained = 'Unstained.fcs', fsc = 'FSC-A', ssc = 'SSC-A',
                         patt = '-A$', method = 'median')   # flowStats: returns the matrix directly
```

## Logicle / Biexponential Transform (fluorescence)

**Goal:** Display and analyze compensated fluorescence with negatives handled honestly.

**Approach:** `estimateLogicle()` (flowWorkspace) derives `w` from the data's most-negative events; apply with `transform()`.

```r
library(flowWorkspace)

fluo <- colnames(fcs_comp)[grepl('-A$', colnames(fcs_comp)) & !grepl('FSC|SSC', colnames(fcs_comp))]
lgcl <- estimateLogicle(fcs_comp, channels = fluo)   # data-driven w; t=262144, m=4.5, a=0 defaults
fcs_t <- transform(fcs_comp, lgcl)
```

## Arcsinh Transform (CyTOF cofactor 5; fluorescence/spectral ~150 or per-channel)

**Goal:** Variance-stabilize mass-cytometry counts (or any pipeline feeding clustering).

**Approach:** `asinh(x/cofactor)`; flowCore's `arcsinhTransform` is `asinh(a + b*x) + c`, so set `b=1/cofactor`. CATALYST `prepData` defaults cofactor=5.

```r
COFACTOR <- 5      # standard CyTOF cofactor, codified in the CATALYST workflow (Nowicka 2017); ~150 for fluorescence

asinhT <- arcsinhTransform(transformationId = 'asinh', a = 0, b = 1/COFACTOR, c = 0)
fcs_t  <- transform(fcs, transformList(marker_channels, asinhT))

# CATALYST path (CyTOF): cofactor=5 default; OVERRIDE for fluorescence/spectral
sce <- CATALYST::prepData(fs, panel, md, transform = TRUE, cofactor = COFACTOR)
```

## Per-Method Failure Modes

### Over-compensation (negative pull-down)
**Trigger:** matrix slope over-estimated from dim controls. **Mechanism:** subtraction overshoots. **Symptom:** negative population pulled below zero, "comma" shape. **Fix:** controls at least as bright as the sample; AutoSpill regression; never hand-tune to flatten spread.

### Wrong logicle width clips negatives
**Trigger:** fixed `w` instead of `estimateLogicle`. **Mechanism:** linear region too narrow. **Symptom:** negative population piled on the axis. **Fix:** estimate `w` on compensated data.

### Cofactor compresses a dim marker
**Trigger:** cofactor 5 on fluorescence (or 150 on CyTOF). **Mechanism:** linear region mismatched to the noise band. **Symptom:** dim-positive collapses into the negative; clusters don't reproduce. **Fix:** 5 for CyTOF, ~150 for fluorescence; per-channel via `flowVS::estParamFlowVS`.

### Compensating spectral data
**Trigger:** treating Aurora data as conventional. **Mechanism:** subtraction is the wrong model for an overdetermined system. **Symptom:** residual spread, false positives. **Fix:** unmix against single-stain reference spectra + unstained AF.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| arcsinh cofactor = 5 (mass) | Nowicka 2017 *F1000Res* 6:748 (CATALYST workflow) | matches CyTOF ion-count near-zero noise band |
| arcsinh cofactor ~150 (fluorescence) | community/CATALYST convention (not a derived optimum) | PMT photon scale is far larger; per-channel flowVS supersedes |
| comp control >= sample brightness | Roederer 2001 *Cytometry* 45:194 | slope estimated over the widest lever arm; extrapolation amplifies error |
| spreading is intensity-dependent (~sqrt of signal) | Nguyen 2013 *Cytometry A* 83:306 | SSM is normalized to be gain-independent for panel design |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `compensate()` channel mismatch | matrix colnames != FCS channels | align names before compensate |
| all-negative after transform | transform applied before/without compensation | compensate on linear data first |
| `estimateLogicle` not found | called from flowCore | it lives in flowWorkspace |
| `arcsinhTransform` ignores "cofactor" | param is `b`, not `cofactor` | set `b = 1/cofactor` |

## References

- Roederer 2001 *Cytometry* 45(3):194-205 — spreading error / compensation artifacts.
- Nguyen 2013 *Cytometry A* 83(3):306-315 — spillover spreading matrix; panel design.
- Roca 2021 *Nat Commun* 12:2890 — AutoSpill robust-regression compensation.
- Parks 2006 *Cytometry A* 69(6):541-551 — logicle display.
- Moore & Parks 2012 *Cytometry A* 81(4):273-277 — logicle operational update.
- Bendall 2011 *Science* 332(6030):687-696 — CyTOF mass cytometry; arcsinh-median analysis.
- Nowicka 2017 *F1000Research* 6:748 — CATALYST workflow; codifies the cofactor-5 convention.
- Azad 2016 *BMC Bioinformatics* 17:291 — flowVS per-channel cofactor.
- Chevrier 2018 *Cell Syst* 6(5):612-620 — CyTOF spillover compensation (CATALYST).

## Related Skills

- fcs-handling - Load FCS and retrieve the spillover keyword first
- gating-analysis - Gate on the transformed, compensated scale
- clustering-phenotyping - Cluster compensated, transformed data
- cytometry-qc - QC before and after preprocessing
- imaging-mass-cytometry/data-preprocessing - Shared arcsinh/metal-channel conventions
- spatial-transcriptomics/spatial-proteomics - Spectral/metal unmixing context

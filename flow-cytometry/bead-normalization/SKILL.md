---
name: bio-flow-cytometry-bead-normalization
description: Bead-based signal normalization and cross-batch harmonization for CyTOF and high-parameter cytometry - EQ four-element bead normalization of instrument sensitivity drift (CATALYST normCytof, premessa), and reference-anchor cross-batch normalization (CytoNorm, per-cluster quantile splines). Covers the distinction between within-run drift correction and between-batch correction, the mandatory anchor/reference sample, why normalization is per-cluster with many quantiles, and the over-correction risk. Use when correcting CyTOF signal drift, harmonizing multi-batch or multi-site studies, or deciding whether to normalize data versus model batch in the design.
tool_type: r
primary_tool: CATALYST
---

## Version Compatibility

Reference examples tested with: CATALYST 1.26+, CytoNorm 2.0+, flowCore 2.14+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

`normCytof()` returns a LIST (`$data`, `$beads`, `$removed`, ...), not a flowFrame; `beads="dvs"` encodes EQ masses 140,151,153,165,175. Confirm with `?normCytof` before relying on slot names.

# Bead Normalization

**"Normalize my CyTOF data"** -> Correct instrument sensitivity drift with EQ beads (within/across runs), then harmonize batches with a reference anchor.
- R (drift): `CATALYST::normCytof()` (EQ-bead-based) or `premessa`
- R (batch): `CytoNorm::CytoNorm.train()` + `CytoNorm.normalize()` (per-cluster quantile splines)

## The Single Most Important Modern Insight -- Two Different Layers; Anchor Controls Are the Guarantee

Bead normalization and batch normalization correct DIFFERENT things and are NOT interchangeable. (1) EQ-BEAD normalization (Finck 2013 *Cytometry A* 83:483) corrects within-run and run-to-run instrument SENSITIVITY DRIFT using the four-element beads as a physical internal standard - applied first, on raw counts. (2) CROSS-BATCH normalization (CytoNorm, Van Gassen 2020 *Cytometry A* 97:268) corrects staining/acquisition batch effects using a shared ANCHOR/reference sample present in EVERY batch, learning per-FlowSOM-cluster quantile-spline transforms. Beads cannot fix staining-batch or reagent-lot effects; CytoNorm cannot fix intra-run detector drift. The anchor control is the load-bearing design element: because it is biologically identical across batches, any cross-batch difference in it is technical BY CONSTRUCTION. Dropping the anchor (CytoNorm 2.0) is convenient but reintroduces the over-correction risk the anchor was designed to eliminate - so the safest stance for inference is to MODEL batch in the diffcyt design and reserve normalization for visualization/clustering display.

## Why Per-Cluster and Many (99) Quantiles

Batch effects are cell-type-specific - a marker can drift in monocytes but not in T cells - so a single global channel transform over-corrects one population while under-correcting another and can erase real abundance differences. CytoNorm therefore learns the transform PER FlowSOM cluster. And it uses ~99 quantiles + a spline because the drift is non-linear and intensity-dependent (the negative and positive peaks move by different amounts); a single median shift or linear rescale reintroduces the distortion it is trying to remove.

## EQ-Bead Normalization (drift)

**Goal:** Correct sensitivity drift and remove bead events.

**Approach:** `normCytof()` gates beads, computes the correction on the linear scale, and returns a list - the cleaned SCE is in `$data`.

```r
library(CATALYST)

sce <- prepData(fs, panel, md)                         # no by_time arg - normalization is normCytof's job
res <- normCytof(sce, beads = 'dvs',                   # EQ masses 140,151,153,165,175
                 k = 500, remove_beads = TRUE, overwrite = FALSE)   # k = smoothing window (default; affects bead-trace viz, not correction magnitude)
sce_norm <- res$data                                   # normalized SCE; res$beads / res$removed available
```

## Cross-Batch Normalization (CytoNorm)

**Goal:** Harmonize batches using a shared reference sample.

**Approach:** Train on the anchor (present in every batch) -> learn per-cluster quantile splines -> apply to the real samples. `testCV()` first: if cluster CV is high, the FlowSOM model is batch-unstable and per-cluster splines will distort (fall back to `nClus=1`).

```r
library(CytoNorm)

model <- CytoNorm.train(files = ref_files, labels = batch_labels, channels = marker_channels,
                        transformList = tl,
                        FlowSOM.params = list(nCells = 6000, xdim = 10, ydim = 10, nClus = 10),
                        normMethod.train = QuantileNorm.train,
                        normParams = list(nQ = 99), seed = 42)
CytoNorm.normalize(model = model, files = sample_files, labels = batch_labels,
                   transformList = tl, transformList.reverse = tl_rev,   # BOTH required
                   outputDir = 'normalized/')
```

## Per-Method Failure Modes

### Treating bead and batch normalization as the same
**Trigger:** expecting beads to fix staining-batch effects. **Mechanism:** different layers. **Symptom:** residual batch structure after bead norm. **Fix:** bead norm for drift; CytoNorm for batch.

### No anchor in a batch
**Trigger:** a batch lacking the reference sample. **Mechanism:** nothing biologically-identical to learn from. **Symptom:** that batch can't be normalized / is over-corrected. **Fix:** run the anchor in every batch (or model batch instead).

### Over-correction
**Trigger:** CytoNorm with groups confounded with batch, or anchor-free on variable samples. **Mechanism:** splines absorb real biology. **Symptom:** attenuated group differences. **Fix:** `testCV()` check; model batch in diffcyt for inference; normalize for display only.

### Using normCytof return as a flowFrame
**Trigger:** `sce_norm <- normCytof(...)`. **Mechanism:** it returns a list. **Symptom:** downstream type error. **Fix:** `res$data`.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| bead drift reduced ~4.9x -> 1.3x | Finck 2013 *Cytometry A* 83:483 | EQ-bead correction over a month of runs |
| 99 quantiles, per-cluster | Van Gassen 2020 *Cytometry A* 97:268 | non-linear intensity-dependent, cell-type-specific drift |
| EQ masses 140,151,153,165,175 (`dvs`) | CATALYST | DVS/Fluidigm EQ four-element bead set |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `normCytof` output not usable | it returns a list | use `res$data` |
| `prepData(by_time=TRUE)` errors | no such argument | use `normCytof()` for bead/drift correction |
| CytoNorm distorts populations | unstable FlowSOM clustering | run `testCV()`; reduce `nClus` (or 1) |
| batch effect remains | only bead-normalized | add CytoNorm with anchor samples |

## References

- Finck 2013 *Cytometry A* 83(5):483-494 — EQ-bead normalization of CyTOF drift.
- Van Gassen 2020 *Cytometry A* 97(3):268-278 — CytoNorm per-cluster quantile normalization.
- Quintelier 2025 *Cytometry A* 107(2):69-87 — CytoNorm 2.0 (anchor-free; over-correction caveat).
- Chevrier 2018 *Cell Syst* 6(5):612-620 — CyTOF spillover (CATALYST normalization context).

## Related Skills

Workflow order (CyTOF): EQ-bead drift normalization (raw counts, FIRST) -> cytometry-qc -> doublet-detection -> clustering -> CytoNorm cross-batch (LAST). The two normalization layers sit at opposite ends.

- cytometry-qc - EQ-bead-median-vs-Time is the primary CyTOF drift readout
- doublet-detection - Remove doublets before normalization
- compensation-transformation - Transform scale used by CytoNorm
- clustering-phenotyping - Cluster across normalized batches
- differential-analysis - Model batch in the design rather than over-cleaning
- experimental-design/batch-design - Anchor/reference-sample design; differential-expression/batch-correction for execution

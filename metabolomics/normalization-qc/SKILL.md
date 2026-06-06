---
name: bio-metabolomics-normalization-qc
description: Designs QC, corrects signal drift, removes batch effects, filters features, normalizes samples, and imputes missing values for untargeted LC-MS/GC-MS metabolomics, framing each step as a measurement model that can create or erase biological signal. Use when processing a peak/feature table before statistical analysis, choosing a drift-correction or sample-normalization method, deciding QC RSD vs D-ratio filtering, or handling left-censored missing values. The feature table is produced by metabolomics/xcms-preprocessing or metabolomics/msdial-preprocessing; transformation/scaling for modeling defers to metabolomics/statistical-analysis; cross-study design issues link to experimental-design/batch-design.
tool_type: r
primary_tool: pmp
---

## Version Compatibility

Reference examples tested with: pmp 1.14+, statTarget 1.30+, imputeLCMD 2.1+, missForest 1.5+, sva 3.50+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

Valid drift correction requires QC injections that bracket the samples at both ends and sample the drift curve (~1 QC every 5-10 injections); conditioning injections must be excluded. Valid batch correction requires biological groups randomized across batches; a confounded design cannot be rescued by any algorithm.

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Metabolomics Normalization and QC

**"Normalize my metabolomics data and correct for batch effects"** -> Filter junk features by QC quality, correct within-batch drift against injection order, normalize per-sample dilution, and impute by missingness mechanism -- each step verified against held-out QCs, not just QC clustering.
- R drift correction: `QCRSC()` (pmp), `shiftCor()` (statTarget)
- R normalization: `pqn_normalisation()` (pmp)
- R imputation: `mv_imputation()` (pmp), `impute.QRILC()` (imputeLCMD), `missForest()` (missForest)

## The Single Most Important Insight -- Normalization Is a Modeling Decision, Not a Cleanup Step

Every preprocessing step imposes an assumption about where the *unwanted* variance lives; if that assumption is wrong the result is not noisier, it is confidently wrong. Three corollaries reorganize the whole skill. (1) QC-based correction assumes the pooled QC's per-feature drift trajectory *is* the samples' trajectory -- false for subgroup-specific features (the pool dilutes them toward absence) and for features at different abundance in samples vs pool (suppression is concentration-dependent), so correcting them extrapolates from noise. (2) The order/batch/biology confound is information-theoretically unwinnable post hoc: if group is collinear with batch or injection order, no estimator can attribute the shared variance to one source -- it only redistributes it, wrongly. Randomization at the bench is the only real fix. (3) Over-correction is invisible to the metric everyone reports: "QC RSD dropped / QCs cluster tighter" is exactly what a too-flexible model games (a cubic spline threading every QC drives QC RSD to ~0% while raising biological-sample RSD). Validate on held-out QCs and dilution-QC linearity, never on the metric the model optimized.

## The Four Orthogonal Operations (Do Not Conflate)

| Operation | Acts on | Removes | Methods |
|---|---|---|---|
| Drift / signal correction | each feature, within a batch, vs injection order | longitudinal intensity decay/rise (column fouling, sensitivity loss) | QC-RLSC (LOESS), QCRSC (spline), QC-RFSC (RF vs order), SERRF (RF across correlated features) |
| Batch correction | each feature, across batches | step-changes between analytical batches | QC-anchored median/reference alignment; ComBat (reserved, dangerous) |
| Sample normalization | each sample (column) | dilution / total-amount differences | PQN, MSTUS, TIC/sum, median, internal standard |
| Transformation / scaling | each feature (row) | mean-variance dependence; range dominance | log/glog; Pareto/auto -> defers to metabolomics/statistical-analysis |

TIC normalization does not handle drift, and -- because of closure -- can spread one feature's change across all others. Keep the axes separate.

## Pipeline Order (and Why Order Matters)

| # | Step | Why here | Tool |
|---|---|---|---|
| 0 | Exclude conditioning injections | Pre-equilibrium signal warps a LOESS edge and corrupts RSD/blank filters | manual (drop first ~8 QC) |
| 1 | Blank filter -> detection-rate filter | Removes background/contaminant and mostly-absent features before any model trains on them | `filter_peaks_by_blank`, `filter_peaks_by_fraction` (pmp) |
| 2 | Within-batch drift correction | Flattens order-dependent trend per feature before cross-sample comparison | `QCRSC` (pmp), `shiftCor` (statTarget) |
| 3 | QC RSD / D-ratio filter | Drift correction *should* improve RSD; filter after so reproducibility reflects corrected data (report both stages) | `filter_peaks_by_rsd` (pmp), `dratio_filter` (structToolbox) |
| 4 | Between-batch alignment | QC-anchored offsets removed after within-batch drift is flat | median-of-QC / batchCorr |
| 5 | Missing-value imputation | Filter aggressively first, then impute only the sparse residual holes by mechanism | `mv_imputation` (pmp), `impute.QRILC`, `missForest` |
| 6 | Sample normalization | Dilution correction on quality features, after junk removed | `pqn_normalisation` (pmp) |
| 7 | Transformation + scaling | Defers to metabolomics/statistical-analysis | `glog_transformation` (pmp) |

Detection-rate filtering must precede imputation: never impute a feature that is 90% missing, which would fabricate 90% of it.

## Decision Tree -- Sample Normalization by Matrix

| Matrix / situation | Use | Why |
|---|---|---|
| Urine / variable-dilution biofluid | PQN or MSTUS (osmolality/SG if measured) | Dilution varies wildly; PQN's median-quotient isolates the common dilution factor; MSTUS excludes drug/diet xenobiotics that corrupt TIC |
| Plasma / serum | PQN or median (TIC only if no dominant peak) | Volume relatively constant; closure risk lower but still present |
| Tissue / cells | Per measured amount (mass, protein, cell count) at the bench | The confounder (input amount) is known -- more honest than any data-driven post-hoc method |
| Targeted / few analytes | Per-class internal standards | One IS cannot represent all chemical classes/RT regions |
| Global profile genuinely differs between groups | Avoid quantile normalization | It forces all samples to one distribution, erasing real distributional biology |
| Creatinine for urine | Avoid as sole method | Fails under renal impairment / muscle-mass differences (Warrack 2009) |

When >50% of features move coherently (potent drug, gross pathology), the PQN median-quotient measures the biology, not dilution, and subtracts it out -- switch to a measured external quantity and check whether the normalization factor correlates with the phenotype.

## Decision Tree -- Drift Correction Method

| Situation | Do | Why |
|---|---|---|
| Smooth monotonic drift, frequent QCs, small/medium study | QCRSC (spline) or QC-RLSC | Per-feature fit vs order; CV-select span to avoid overfit |
| Non-smooth / multi-pattern drift within a batch | QC-RFSC (statTarget) or batchCorr clusters | RF / cluster-based captures non-monotonic trend |
| Large cohort (>~500), complex multi-source error, want lowest RSD | SERRF | Borrows strength across correlated features (~5% RSD on >800-sample cohorts, Fan 2019) |
| Sparse QCs (<5-6 spanning the batch) | Coarse median-of-QC offset or no within-batch correction | LOESS/spline with too few QCs produces gaps/garbage |
| Feature weak/absent in QCs | Exclude from correction | Correcting it extrapolates from noise |
| No detectable drift in a feature | Do not correct it | Correcting a flat QC trajectory only adds the model's wiggle |
| Run order confounded with biology | Do not drift-correct; fix design or caveat | A smooth function of order absorbs and subtracts the biological trend |

Flexible ML methods (SERRF/RF/adversarial) win on large complex cohorts but are *more* prone to learning-and-removing biology that tracks order/batch. Always confirm QC RSD dropped AND biological-sample RSD did not rise.

## Filter Features by QC Quality (RSD and D-ratio)

**Goal:** Keep only reproducible features whose technical variance is small relative to biological variance.

**Approach:** Compute per-feature QC RSD and the robust D-ratio (technical SD / biological SD), then apply a boolean mask. Lead with D-ratio: CV alone is matrix-blind, scoring a precisely-measured-but-flat feature as good and a noisy-but-biologically-huge feature as bad.

```r
library(matrixStats)

robust_dratio_filter <- function(data, is_qc, dratio_max = 0.5, rsd_max = 0.3) {
    qc <- as.matrix(data[is_qc, ])
    bio <- as.matrix(data[!is_qc, ])
    # MAD-based (robust) form, because MS intensities are right-skewed
    sd_qc <- colMads(qc, na.rm = TRUE)
    sd_bio <- colMads(bio, na.rm = TRUE)
    dratio <- sd_qc / sd_bio
    rsd <- colSds(qc, na.rm = TRUE) / colMeans(qc, na.rm = TRUE)
    keep <- dratio <= dratio_max & rsd <= rsd_max
    keep[is.na(keep)] <- FALSE
    message(sprintf('D-ratio<=%.2f & RSD<=%.0f%%: kept %d / %d features',
                    dratio_max, rsd_max * 100, sum(keep), ncol(data)))
    data[, keep]
}
```

## Correct Within-Batch Drift (QC-RSC)

**Goal:** Flatten per-feature, injection-order-dependent signal drift using the QC trajectory.

**Approach:** Fit a QC-robust smoothing spline of intensity vs injection order per feature, interpolate at every sample position, and divide. pmp's `QCRSC` selects the spline smoothing by leave-one-out CV when `spar=0`, requires `minQC` QCs per batch, and excludes features too weak in QC automatically.

```r
library(pmp)

# df: features in ROWS, samples in COLUMNS (pmp convention)
corrected <- QCRSC(df = feature_matrix, order = injection_order, batch = batch_id,
                   classes = sample_class, spar = 0, log = TRUE,
                   minQC = 5, qc_label = 'QC')
# Verify correction worked on HELD-OUT QCs / dilution linearity, not on QC clustering.
```

statTarget alternative (`MLmethod='QCRFSC'` for RF, `'QCRLSC'` for LOESS; `QCspan=0` auto-GCV span applies to QCRLSC; inputs are two order-aligned CSVs):

```r
library(statTarget)
shiftCor(samPeno = 'meta.csv', samFile = 'peaks.csv', Frule = 0.8,
         MLmethod = 'QCRFSC', ntree = 500, QCspan = 0, degree = 2,
         imputeM = 'KNN', coCV = 30, plot = FALSE)
```

## Normalize Per-Sample Dilution (PQN)

**Goal:** Remove per-sample global intensity differences (dilution, extraction efficiency) without subtracting genuine fold changes.

**Approach:** Build a reference spectrum (median of QCs), compute per-feature sample/reference quotients, take the median quotient as the dilution factor, and divide. The median is robust because it ignores the minority of genuinely-changed features.

```r
library(pmp)
# df: features in ROWS, samples in COLUMNS; reference built from QC samples
normalized <- pqn_normalisation(df = feature_matrix, classes = sample_class,
                                qc_label = 'QC')
```

## Impute by Missingness Mechanism

**Goal:** Fill residual sparse holes with the method matched to *why* the value is missing.

**Approach:** Diagnose the mechanism per feature -- missingness correlated with low abundance is MNAR (left-censored) and needs QRILC/GSimp; sporadic missingness across the abundance range is MAR and needs RF/kNN. Using a MAR method on MNAR zeros pulls the censored group's mean up and erases the on/off signal.

```r
library(imputeLCMD)
library(missForest)

# MNAR / left-censored: random draws from a fitted truncated-normal (features in ROWS)
qrilc_imputed <- impute.QRILC(feature_matrix_features_in_rows, tune.sigma = 1)[[1]]

# MAR / sporadic: iterative random-forest prediction (samples in ROWS, features in COLS)
rf_imputed <- missForest(sample_by_feature_matrix, maxiter = 10, ntree = 100)$ximp
```

Half-min imputation collapses the imputed subset's variance to zero, understating SE and inflating false significance -- prefer QRILC/GSimp, which draw a distribution of plausible low values. Re-run key results under >=2 imputation methods; if headline metabolites flip, the finding lives in the imputation.

## Per-Method Failure Modes

### QC-spline / LOESS over-correction
- **Trigger:** Span too small, too few QCs, or a cubic spline using each QC as a lock-point.
- **Mechanism:** The fit threads the QCs perfectly, modeling inter-QC noise as drift and injecting it into samples.
- **Symptom:** QC RSD drops toward 0% while biological-sample RSD rises and the number of significant features explodes.
- **Fix:** CV-select the span, validate on held-out QCs and dilution-QC linearity, compare biological-sample RSD before/after, not just QC RSD.

### Edge extrapolation
- **Trigger:** Samples injected before the first QC or after the last QC.
- **Mechanism:** LOESS/spline interpolate between QCs but extrapolate beyond them, producing unstable correction factors.
- **Fix:** Bracket samples with QCs at both ends; exclude conditioning injections.

### TIC closure artifact
- **Trigger:** One dominant or up-regulated feature under TIC/sum normalization.
- **Mechanism:** The constant-sum constraint forces every other feature's normalized value down, manufacturing apparent coordinated down-regulation (Aitchison closure; spurious negative correlations).
- **Fix:** Use PQN/MSTUS or log-ratios; if "many features moved together," suspect closure from one big mover before believing coordination.

### ComBat under imbalance
- **Trigger:** Biological groups unbalanced across batches.
- **Mechanism:** Empirical-Bayes shrinkage confounds class with batch; under imbalance it can fabricate thousands of false differences or, without the covariate, delete real ones (Nygaard 2016).
- **Fix:** Prefer QC-anchored between-batch alignment (the pool has no group, so it cannot confound). Reserve ComBat for balanced designs and always pass the biological covariate via `mod=`. Randomize so it is never needed.

### Wrong imputation mechanism
- **Trigger:** kNN/RF applied to below-LOD (MNAR) values.
- **Mechanism:** MAR methods borrow abundance from detected samples, pulling the censored group's mean up and shrinking the very difference under test.
- **Fix:** Diagnose mechanism per feature; QRILC/GSimp for left-censored, RF/kNN only for sporadic MAR.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|---|---|---|
| QC RSD <= 20-30% (15% gold) | Dunn 2011; Broadhurst 2018 | Reproducibility floor in the matrix-matched pool; 20% aspirational for LC-MS, 30% common |
| D-ratio <= 0.5 (0.2 excellent), robust/MAD form | Broadhurst 2018 | Technical SD < biological SD -- the honest, matrix-aware filter; MAD form because MS intensities are right-skewed |
| Blank ratio >= 3-5x | community convention (Dunn lineage) | Features below 3-5x blank are dominated by background/carryover, not biology |
| Detection rate >= 50-80% (or 80% within any one group) | statTarget Frule=0.8; pmp filter_peaks_by_fraction | Reliable signal; "within any group" preserves on/off group-specific metabolites |
| Dilution-QC correlation r >= 0.7-0.8 | community convention | Real metabolites scale with dilution; artefacts/in-source ions do not |
| QCs >= 5-10% of injections, ~1 every 5-10 samples | Broadhurst 2018; mQACC 2022 | Must sample the drift curve densely enough to avoid LOESS extrapolation |

Thresholds are conventions, not laws: choose them a priori, report each one, and report how many features each filter removed (mQACC reporting standard).

## Common Errors

| Error / symptom | Cause | Solution |
|---|---|---|
| `could not find function "statTarget"` | No such entry point | Use `shiftCor()` (drift correction) and `statAnalysis()` (post-hoc stats) |
| `mv_imputation` errors on `method='sm'` | Small-value method is `'sv'`, not `'sm'` | Use `method='sv'` (also valid: `knn`, `rf`, `bpca`, `mn`, `md`) |
| QRILC output is malformed | `impute.QRILC` returns a list, not a matrix; expects features in rows | Index `[[1]]`; transpose so features are rows |
| MetaboAnalystR `Normalization` errors | `SanityCheckData(mSet)` not run first | Call `SanityCheckData` -> `ReplaceMin` -> `Normalization` in order |
| Correction made data worse | Span overfit / weak-in-QC features corrected / order confounded with biology | Back off span, exclude weak-in-QC features, check randomization |
| Effect vanished after drift correction | Run order confounded with group; trend absorbed the biology | Check the design; report drift and effect as inseparable if confounded |

## References

- Dunn WB, Broadhurst D, Begley P, et al. 2011. Procedures for large-scale metabolic profiling of serum and plasma using gas chromatography and liquid chromatography coupled to mass spectrometry. *Nature Protocols* 6:1060-1083.
- Broadhurst D, Goodacre R, Reinke SN, et al. 2018. Guidelines and considerations for the use of system suitability and quality control samples in mass spectrometry assays applied in untargeted clinical metabolomic studies. *Metabolomics* 14:72.
- Dieterle F, Ross A, Schlotterbeck G, Senn H. 2006. Probabilistic Quotient Normalization as Robust Method to Account for Dilution of Complex Biological Mixtures. Application in 1H NMR Metabonomics. *Analytical Chemistry* 78:4281-4290.
- Fan S, Kind T, Cajka T, et al. 2019. Systematic Error Removal Using Random Forest for Normalizing Large-Scale Untargeted Lipidomics Data. *Analytical Chemistry* 91:3590-3596.
- Brunius C, Shi L, Landberg R. 2016. Large-scale untargeted LC-MS metabolomics data correction using between-batch feature alignment and cluster-based within-batch signal intensity drift correction. *Metabolomics* 12:173.
- Luan H, Ji F, Chen Y, Cai Z. 2018. statTarget: A streamlined tool for signal drift correction and interpretations of quantitative mass spectrometry-based omics data. *Analytica Chimica Acta* 1036:66-72.
- Wei R, Wang J, Su M, et al. 2018. Missing Value Imputation Approach for Mass Spectrometry-based Metabolomics Data. *Scientific Reports* 8:663.
- Wei R, Wang J, Jia E, et al. 2018. GSimp: A Gibbs sampler based left-censored missing value imputation approach for metabolomics studies. *PLoS Computational Biology* 14:e1005973.
- Stekhoven DJ, Buhlmann P. 2012. MissForest -- non-parametric missing value imputation for mixed-type data. *Bioinformatics* 28:112-118.
- Nygaard V, Rodland EA, Hovig E. 2016. Methods that remove batch effects while retaining group differences may lead to exaggerated confidence in downstream analyses. *Biostatistics* 17:29-39.
- van den Berg RA, Hoefsloot HCJ, Westerhuis JA, et al. 2006. Centering, scaling, and transformations: improving the biological information content of metabolomics data. *BMC Genomics* 7:142.
- Warrack BM, Hnatyshyn S, Ott KH, et al. 2009. Normalization strategies for metabonomic analysis of urine samples. *Journal of Chromatography B* 877:547-552.
- Thonusin C, IglayReger HB, Soni T, et al. 2017. Evaluation of intensity drift correction strategies using MetaboDrift, a normalization tool for multi-batch metabolomics data. *Journal of Chromatography A* 1523:265-274.
- Chamberlain CA, Rubio VY, Garrett TJ. 2019. Impact of matrix effects and ionization efficiency in non-quantitative untargeted metabolomics. *Metabolomics* 15:135.
- Wehrens R, Hageman JA, van Eeuwijk F, et al. 2016. Improved batch correction in untargeted MS-based metabolomics. *Metabolomics* 12:88.

## Related Skills

- metabolomics/xcms-preprocessing - Generates the feature table this skill consumes
- metabolomics/msdial-preprocessing - Alternative feature-table source
- metabolomics/statistical-analysis - Transformation/scaling and downstream multivariate stats
- experimental-design/batch-design - Randomization and design that make correction valid
- differential-expression/batch-correction - ComBat/SVA mechanics shared with transcriptomics

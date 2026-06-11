---
name: bio-multi-omics-mixomics-analysis
description: Builds supervised and unsupervised multivariate integration across bulk omics blocks with mixOmics - sPLS for sparse pairwise correlation, DIABLO (block.splsda) for a multi-block discriminant signature, rCCA for regularized canonical correlation, and MINT for multi-study integration. Covers why these projection methods maximize covariance or correlation and not truth, why DIABLO's design matrix is the central correlation-versus-discrimination decision, why cross-validation must wrap keepX selection or the reported error is leaked, why balanced error rate is required under class imbalance, and why DIABLO needs matched samples while MINT handles multiple cohorts. Use when finding a cross-omic discriminant signature for a known outcome, selecting correlated features between two omics, tuning keepX, or integrating one omic across studies. For unsupervised factors see mofa-integration; for the method decision see integration-design; for cross-validation theory see machine-learning/model-validation.
tool_type: r
primary_tool: mixOmics
---

## Version Compatibility

Reference examples tested with: mixOmics 6.26+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('mixOmics')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

The DIABLO model is `block.splsda` - there is no function literally named `diablo`. mixOmics input is samples-by-features (the opposite of MOFA2's features-by-samples), and `scale=TRUE` is the default, so feed appropriately transformed but not pre-standardized matrices.

# mixOmics Multi-Omics Analysis

**"Find a cross-omic signature that discriminates my groups"** -> Project the blocks onto sparse latent components that maximize an association criterion - because these methods maximize covariance, not truth, so a signature that discriminates the training cohort is guaranteed and only out-of-sample replication makes it real.
- R: `block.splsda()` (DIABLO, supervised), `spls()` (pairwise), `mint.splsda()` (multi-study)

Scope: supervised multi-block discriminant integration (DIABLO), sparse pairwise correlation (sPLS), regularized CCA (rCCA), and multi-study integration (MINT). Unsupervised factor models -> mofa-integration. The method-selection decision -> integration-design. Generic cross-validation / leakage theory -> machine-learning/model-validation. Per-omic scaling/batch -> data-harmonization. Enrichment of selected features -> pathway-analysis/go-enrichment.

## The Single Most Important Modern Insight -- DIABLO Maximizes Covariance, Not Truth, So a Signature Scored on the Data That Selected It Is an Artifact

These are projection methods: they find linear combinations of features that optimize an association criterion (covariance for PLS/DIABLO, correlation for CCA), and DIABLO additionally optimizes whatever the design matrix tells it to. With n much smaller than p the methods can fit almost anything, so the output is never "the cross-omic drivers of the disease" - it is the features that best satisfied the chosen criterion on the available samples. Three rules follow:

1. **Cross-validation must wrap feature selection, not follow it.** Tuning keepX on the full data and then reporting `perf()` cross-validation error on the same data is leakage (up to ~0.15 AUC inflation). mixOmics selects inside folds when tuning, which is correct for choosing keepX - but the minimized CV error it returns is still optimistic. The honest number comes from an external test set never touched during tuning, or a fully nested CV.
2. **The design matrix is the central decision, not a default.** The off-diagonal weights in [0,1] trade discrimination against cross-block correlation: near 1 gives biologically coherent, inter-correlated signatures at the cost of classification; near 0 gives the best classification with a disconnected network. The tutorials' 0.1 leans toward prediction and is not a recommendation from the DIABLO paper - choose it from the goal and report it.
3. **The selected features are candidates, not biomarkers.** selectVar discriminates the training cohort by construction; that is the null result, not the finding. Replication in an independent cohort (or MINT across studies) is the finding.

## Tool Taxonomy

| Method (mixOmics fn) | Citation | Optimizes | Supervision |
|----------------------|----------|-----------|-------------|
| sPLS (`spls`) | Le Cao 2008 *Stat Appl Genet Mol Biol* 7:Article 35 | covariance between TWO blocks, sparse | unsupervised pairwise |
| rCCA (`rcc`) | Gonzalez 2008 *J Stat Softw* 23(12) | correlation between two blocks, ridge/shrinkage regularized | unsupervised pairwise |
| DIABLO (`block.splsda`) | Singh 2019 *Bioinformatics* 35:3055 | covariance across MANY blocks + discrimination, sparse | supervised, multi-block |
| MINT (`mint.splsda`) | Rohart 2017 *BMC Bioinformatics* 18:128 | one omic across studies, study as fixed effect | supervised, horizontal |
| sPLS-DA (`splsda`) | Le Cao 2011 *BMC Bioinformatics* 12:253 | discrimination in ONE block, sparse | supervised, single-block |
| sPCA (`spca`) | Shen 2008 *J Multivar Anal* 99:1015 | variance in one block, sparse | unsupervised, single-block |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Two or more omics, matched samples, a categorical outcome | DIABLO (`block.splsda`) | supervised multi-block; design tunes correlation vs discrimination |
| Two omics, no outcome, want a sparse correlated feature list | sPLS canonical (`spls`, mode='canonical') | covariance, symmetric, feature selection |
| Two omics, no outcome, want the global correlation landscape | rCCA (`rcc`, shrinkage) | correlation criterion, regularized for p>n |
| One omic predicts another (directional) | sPLS regression (`spls`, mode='regression') | asymmetric, Y as response |
| SAME omic across multiple cohorts, reproducible signature | MINT (`mint.splsda`) | horizontal; study as a fixed effect |
| No outcome, want variance-attributed factors, missing blocks OK | -> mofa-integration | Bayesian factor model, not a projection |
| Need an honest performance number | external test set or nested CV | the tuned CV error is optimistic |
| Which method at all / paired vs horizontal | -> integration-design | the correspondence and supervision decision |

## Set Up Matched Blocks

**Goal:** Guarantee the row-by-row sample matching DIABLO/sPLS/rCCA require, so the model correlates feature vectors of the same individuals.

**Approach:** Intersect to common samples and verify identical rowname order across every block; for combining one omic across cohorts, that is horizontal integration and needs MINT, not DIABLO.

```r
library(mixOmics)

common <- Reduce(intersect, list(rownames(X_rna), rownames(X_prot)))   # samples x features
X_blocks <- list(RNA=X_rna[common, ], Protein=X_prot[common, ])
Y <- factor(pheno[common, 'Condition'])
stopifnot(identical(rownames(X_blocks$RNA), rownames(X_blocks$Protein)))   # matched rownames or the result is garbage
```

## sPLS: Sparse Pairwise Correlation

**Goal:** Select a sparse set of features that covary between two omics.

**Approach:** Choose `mode` deliberately - 'canonical' for two omics on equal footing (the CCA-like symmetric framing), 'regression' (the default) only when one block is a designated response. Tune component count, then fit with keepX/keepY feature selection.

```r
spls_res <- spls(X_blocks$RNA, X_blocks$Protein, ncomp=3, mode='canonical',   # symmetric; default 'regression' treats Protein as a response of RNA
                 keepX=c(50, 50, 50), keepY=c(30, 30, 30))
plotVar(spls_res, comp=c(1, 2))
```

## DIABLO: Supervised Multi-Block Signature

**Goal:** Find a cross-omic feature signature that discriminates a known outcome, with the correlation-versus-discrimination trade-off chosen and reported.

**Approach:** Set the design matrix from the goal (high off-diagonal for coherent networks, low for prediction), tune the component count on a non-sparse model, tune keepX inside cross-validation folds using balanced error rate, fit, then report performance on an external set.

```r
design <- matrix(0.5, nrow=length(X_blocks), ncol=length(X_blocks),
                 dimnames=list(names(X_blocks), names(X_blocks)))   # 0.5-1 favors cross-block correlation; <0.5 favors prediction - choose from the goal
diag(design) <- 0

ncomp_fit <- perf(block.plsda(X_blocks, Y, ncomp=5, design=design),               # tune ncomp on a NON-sparse model first
                  validation='Mfold', folds=10, nrepeat=10)                       # nrepeat>=10; a single split is noise

tune <- tune.block.splsda(X_blocks, Y, ncomp=2, design=design,
                          test.keepX=list(RNA=c(10, 25, 50), Protein=c(10, 25, 50)),
                          validation='Mfold', folds=10, nrepeat=10,
                          measure='BER', BPPARAM=BiocParallel::MulticoreParam(workers=4))   # BER, not overall error, for imbalanced classes; cpus= is defunct, use BPPARAM
diablo <- block.splsda(X_blocks, Y, ncomp=2, keepX=tune$choice.keepX, design=design)
```

The features chosen here discriminate the training cohort by construction. For an honest accuracy, hold out an external test set never used in tuning and report `predict()` / `auroc()` on it; the `perf()` error on the tuning data is optimistic because keepX was chosen to minimize it.

## Interpret and Validate

**Goal:** Extract the signature as candidates and visualize cross-block structure without overclaiming.

**Approach:** Pull selected variables per block and component, inspect inter-block correlations, and frame the list as cohort-specific candidates requiring replication.

```r
sel_rna  <- selectVar(diablo, block='RNA', comp=1)$RNA$name        # candidate features, not validated biomarkers
circosPlot(diablo, cutoff=0.7)                                     # inter-block correlations of the selected features
auc <- auroc(diablo, roc.block='RNA', roc.comp=1)
```

## MINT: Multi-Study (Horizontal) Integration

**Goal:** Build a signature for ONE omic that replicates across cohorts by modeling study as a known effect.

**Approach:** Pass a study factor so the model accounts for study-specific variation; this is horizontal integration (same features, different cohorts), the opposite of DIABLO's vertical matched-sample design.

```r
mint_res <- mint.splsda(X=X_rna, Y=Y, study=study, ncomp=3, keepX=c(50, 50, 50))   # study = fixed effect; one omic, many cohorts
plotIndiv(mint_res, study='global', legend=TRUE)
```

## Per-Method Failure Modes

### CV that follows selection instead of wrapping it
**Trigger:** tuning keepX on all samples, then reporting `perf()` CV error on all samples. **Mechanism:** the features were chosen with knowledge of every sample, so no fold is truly held out. **Symptom:** an excellent CV error that collapses in a new cohort. **Fix:** external test set or nested CV; the number used to pick keepX is not an estimate of performance.

### Design matrix copied as a default
**Trigger:** using 0.1 (or any value) without choosing it. **Mechanism:** the off-diagonal trades discrimination against cross-block correlation. **Symptom:** a result marketed as "integrated" while the design told the model to ignore most cross-block correlation. **Fix:** choose the design from the goal, justify it, and ideally show the result under a high and a low weight.

### Un-regularized CCA on p>n
**Trigger:** running plain CCA on omics. **Mechanism:** CCA divides out variances and needs to invert a singular covariance, so it reports correlation 1.0 by overfitting. **Symptom:** perfect, meaningless canonical correlations. **Fix:** use `rcc` with ridge or shrinkage regularization, or sPLS canonical.

### Unmatched samples (or DIABLO where MINT is needed)
**Trigger:** partially overlapping or mis-ordered samples across blocks, or DIABLO on two cohorts of one omic. **Mechanism:** DIABLO/sPLS relate samples row-by-row. **Symptom:** garbage signatures from correlating different individuals. **Fix:** intersect and verify identical rowname order; use MINT for one omic across cohorts.

### Overall error under class imbalance
**Trigger:** tuning on overall classification error with imbalanced classes. **Mechanism:** the majority class dominates the metric. **Symptom:** high accuracy while the minority class is mis-predicted. **Fix:** `measure='BER'` and report per-class error.

### selectVar list presented as validated biomarkers
**Trigger:** calling the selected features mechanistic drivers. **Mechanism:** they discriminate the training cohort by construction. **Symptom:** a biomarker claim that fails to replicate. **Fix:** frame as cohort candidates; require external/cross-study replication.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Design off-diagonal: ~1 for coherence, <0.5 for prediction | Singh 2019 *Bioinformatics* 35:3055 | weights trade cross-block correlation against discrimination; 0.1 is tutorial convention only |
| `nrepeat` >= 10 (50 for a headline number) | mixOmics docs | a single M-fold split is high-variance; nrepeat=1 is illustration only |
| `folds` 5-10 in M-fold CV | mixOmics docs | balances bias and variance of the CV estimate at small n |
| `measure='BER'` for imbalanced classes | mixOmics docs | overall error is dominated by the majority class |
| `ncomp` 1-3, set by `perf()` elbow | Singh 2019 *Bioinformatics* 35:3055 | more components rarely help and risk overfit |
| External test set or nested CV for reported accuracy | machine-learning/model-validation | the tuned CV error is optimistically biased |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `could not find function "diablo"` | DIABLO is `block.splsda` | call `block.splsda`, not `diablo` |
| Perfect canonical correlations | un-regularized CCA on p>n | use `rcc` ridge/shrinkage |
| Reported accuracy fails in a new cohort | CV scored on tuning data | external test set or nested CV |
| `tune.block.splsda` very slow | grid too large / not parallelized | shrink `test.keepX`, set `BPPARAM` (cpus= is defunct) |
| High accuracy, minority class missed | overall error under imbalance | `measure='BER'` |
| Nonsense signature | unmatched/mis-ordered samples | intersect and verify rowname order |

## References

- Le Cao K-A, Rossouw D, Robert-Granie C, Besse P. 2008. A sparse PLS for variable selection when integrating omics data. *Stat Appl Genet Mol Biol* 7:Article 35.
- Le Cao K-A, Boitard S, Besse P. 2011. Sparse PLS discriminant analysis: biologically relevant feature selection and graphical displays for multiclass problems. *BMC Bioinformatics* 12:253.
- Gonzalez I, Dejean S, Martin PGP, Baccini A. 2008. CCA: an R package to extend canonical correlation analysis. *J Stat Softw* 23(12):1-14.
- Rohart F, Gautier B, Singh A, Le Cao K-A. 2017. mixOmics: an R package for 'omics feature selection and multiple data integration. *PLoS Comput Biol* 13:e1005752.
- Rohart F, Eslami A, Matigian N, Bougeard S, Le Cao K-A. 2017. MINT: a multivariate integrative method to identify reproducible molecular signatures across independent experiments and platforms. *BMC Bioinformatics* 18:128.
- Singh A, Shannon CP, Gautier B, et al. 2019. DIABLO: an integrative approach for identifying key molecular drivers from multi-omics assays. *Bioinformatics* 35:3055-3062.
- Shen H, Huang JZ. 2008. Sparse principal component analysis via regularized low rank matrix approximation. *J Multivar Anal* 99:1015-1034.

## Related Skills

- integration-design - The method-selection and paired-vs-horizontal decision
- mofa-integration - Unsupervised factor alternative where no outcome drives the fit
- data-harmonization - Per-block scaling and batch before matrices enter mixOmics
- machine-learning/model-validation - Nested cross-validation and data-leakage theory
- machine-learning/biomarker-discovery - Biomarker-panel selection and validation
- pathway-analysis/go-enrichment - Enrichment of the selected features
- differential-expression/de-results - Single-omic differential expression
- workflows/multi-omics-pipeline - End-to-end multi-omics integration pipeline

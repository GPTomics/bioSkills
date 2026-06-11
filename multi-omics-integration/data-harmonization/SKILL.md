---
name: bio-multi-omics-data-harmonization
description: Harmonizes already-normalized per-omic matrices onto a common footing before joint integration - assembling a MultiAssayExperiment, choosing the per-omic variance-stabilizing transform, deciding per-view versus per-feature scaling, picking a cross-omic batch strategy, and triaging missing data (feature, value, or whole sample; MAR versus MNAR). Covers why a shared-latent integrator is blind to what an omic is so scaling silently decides which block dominates, why batch confounded with biology is irrecoverable and should be modeled as a covariate not scrubbed, and why stacking blocks and running one ComBat erases cross-omic signal. Use when preparing two or more omics for MOFA2, mixOmics, or SNF, deciding a transform or scaling, correcting batch across modalities, or handling missing omics per sample. For deep per-omic normalization see differential-expression, methylation-analysis, proteomics, metabolomics; for the method decision see integration-design; for fusion see mofa-integration, mixomics-analysis.
tool_type: r
primary_tool: MultiAssayExperiment
---

## Version Compatibility

Reference examples tested with: MultiAssayExperiment 1.36+, SummarizedExperiment 1.40+, sva 3.50+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Deep per-omic normalization (DESeq2 VST internals, methylation noob/BMIQ, proteomics VSN) is owned by the per-omic categories and pinned to their own tool versions; this skill calls the cross-omic container and batch tools, so MultiAssayExperiment and sva are the binding versions here.

# Data Harmonization for Multi-Omics

**"Get my omics onto a common footing for integration"** -> Transform each block to a comparable scale, equalize block contribution, decide a single batch strategy, and triage missingness - because the integrator sees one stacked matrix of numbers and spends its first factors on whichever block has the largest raw variance, which is a property of the assay, not the biology.
- R: `MultiAssayExperiment` to coordinate assays, sample map, and colData

Scope: cross-omic harmonization of already-normalized blocks - container assembly, transform choice, per-view scaling, batch strategy, missing-data triage. Deep single-omic normalization -> differential-expression, methylation-analysis, proteomics, metabolomics. The method-selection decision -> integration-design. Fusion math -> mofa-integration, mixomics-analysis. Single-omic RNA batch mechanics -> differential-expression/batch-correction.

## The Single Most Important Modern Insight -- A Shared-Latent Integrator Is Blind to What an Omic Is

A factor model sees one stacked matrix and spends its leading factors on whichever block has the largest raw variance and the most non-Gaussian structure - a property of the measurement scale, not the biology. Harmonization is the act of making the variance budget mean something biological before the model gets to it. Every choice is therefore a silent vote on which omic dominates the shared latent space and which signal is allowed to appear:

1. **The transform decides distribution.** Raw NB counts (0-10^5), bounded methylation betas, right-skewed proteomics, and closed compositional proportions are not comparable; feed them together and the factorization is driven by scale artifacts. Apply the per-omic variance-stabilizing transform FIRST (RNA -> VST/logCPM; methylation -> M-values; proteomics/metabolomics -> log2; compositional -> CLR).
2. **The scaling decides dominance.** Even after transforming, blocks differ in total variance and feature count, and variance is additive across features, so a 20k-gene block out-votes a 200-metabolite block. Per-VIEW scaling equalizes blocks; per-FEATURE unit-variance scaling (mixOmics default) inflates near-constant noise features into spurious factors - filter them first.
3. **The batch model decides what is deleted.** A batch-corrected matrix is the data conditioned on a model of what batch is. When batch is confounded with biology the correction deletes biology with confidence; even on a balanced design, scrub-then-test exaggerates significance (Nygaard 2016). Correct once, prefer modeling batch as a covariate over scrubbing.

The deliverable is not "the harmonized data" - it is a documented chain of transform -> scaling -> batch -> missing-value triage, each justified, because each silently determines the integration result.

## Per-Omic Transform Decision Table

| Omic | Raw-data pathology | Transform for integration | Why |
|------|--------------------|---------------------------|-----|
| RNA-seq (bulk) | NB counts; variance grows with mean; loads on high-count genes | VST or rlog (DESeq2), or log-CPM (voom) | variance-stabilize to approximately homoscedastic and Gaussian (Love 2014); `blind=TRUE` for unsupervised integration |
| Methylation | beta in [0,1]; variance compressed near 0 and 1 | M-value = log2(beta/(1-beta)) for modeling; beta for interpretation | logit un-compresses the extremes -> homoscedastic (Du 2010) |
| Proteomics / metabolomics | positive, right-skewed, multiplicative error; MNAR below detection | log2, then MNAR-aware imputation | log makes multiplicative error additive and Gaussian (Lazar 2016) |
| Microbiome / compositional | simplex, fixed-sum; spurious negative correlation | CLR after zero replacement | maps the simplex to Euclidean space; raw proportions invalid for L2 methods (Gloor 2017) -> metagenomics/abundance-estimation |
| ALL, after the above | blocks differ in total variance and feature count | per-VIEW scaling (MFA singular-value weighting; MOFA `scale_views`) | equalize block contribution without inflating individual features |

The transform is per-omic (owned by the per-omic categories); the scaling is cross-block (owned here). They are sequential and both required - scaling a heteroscedastic block does not make it homoscedastic.

## Harmonization Tool Taxonomy

| Tool | Citation | Role | When |
|------|----------|------|------|
| MultiAssayExperiment | Ramos 2017 *Cancer Res* 77:e39 | container: assays + sampleMap + colData | always - makes sample linkage a structural invariant, not string-munging |
| sva ComBat | Johnson 2007 *Biostatistics* 8:118 | empirical-Bayes batch adjustment (transformed data) | scrub batch PER OMIC when the integrator needs clean input |
| ComBat-seq | Zhang 2020 *NAR Genom Bioinform* 2:lqaa078 | batch adjustment on RNA-seq COUNTS | batch-correct counts before VST, not after |
| limma removeBatchEffect | Ritchie 2015 *Nucleic Acids Res* 43:e47 | regress out batch for visualization | ordination/QC plots only - never feed to an inferential test |
| sva / RUV | Leek 2007 *PLoS Genet* 3:e161; Risso 2014 *Nat Biotechnol* 32:896 | estimate hidden/unwanted variation as covariates | batch is unknown or driven by control genes |
| imputeLCMD / DEP | Lazar 2016 *J Proteome Res* 15:1116 | QRILC/MinProb (MNAR) + kNN (MAR) imputation | below-detection proteomics/metabolomics gaps |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Blocks on different scales going into MOFA/PLS | per-omic transform then per-view scaling | the integrator assumes comparable, homoscedastic inputs |
| One omic has far more features than another | per-view scaling and/or filter the wider view harder | feature count buys variance votes; equalize block contribution |
| Batch differs across omics, integrator needs clean input | ComBat PER OMIC (ComBat-seq for counts) | scrub once per modality; never stack blocks and ComBat together |
| Downstream step is an inferential test | model batch as a covariate (`~ batch + condition`) | scrub-then-test exaggerates confidence (Nygaard 2016) |
| Batch correlates with the condition | do NOT correct; redesign or report the confound | confounded batch is irrecoverable; correcting deletes biology |
| Proteomics missing below detection (MNAR) | QRILC / MinProb (imputeLCMD) | MAR imputation biases low-abundance proteins upward |
| Sporadic within-feature gaps (MAR) | kNN / missForest | local borrowing is valid when missingness is random |
| Some samples missing a whole omic (mosaic) | -> mofa-integration (models the missingness) | do not impute a whole block; intersecting loses scarce n |
| Need the method decision, not the prep | -> integration-design | which integration paradigm fits the question |

## Assemble the Container and Decide Correspondence

**Goal:** Make sample linkage across omics a structural invariant and quantify how mosaic the cohort is, so the impute-versus-model-the-missingness decision is explicit.

**Approach:** Build a MultiAssayExperiment from per-omic SummarizedExperiments; the sampleMap (assay, primary, colname) links assay columns to subjects. `intersectColumns` gives complete cases (subjects), which is what integration wants; `intersectRows` aligns features and is a trap on heterogeneous omics.

```r
library(MultiAssayExperiment)

rna  <- SummarizedExperiment(assays=list(vst=vst_rna), colData=sample_info)       # already VST-normalized
prot <- SummarizedExperiment(assays=list(log2=norm_prot), colData=sample_info)    # already log2 + median-normalized
meth <- SummarizedExperiment(assays=list(mval=m_values), colData=sample_info)     # already M-values

mae <- MultiAssayExperiment(experiments=ExperimentList(RNA=rna, Protein=prot, Methylation=meth),
                            colData=sample_info)
table(complete.cases(mae))        # subjects with every omic
paired <- intersectColumns(mae)   # complete-case fallback; counts the n it costs
```

## Per-View Scaling (Equalize Block Contribution)

**Goal:** Stop the highest-variance or highest-dimensional omic from hijacking the shared factors without inflating noise features.

**Approach:** Filter near-constant features per block, then scale each block so its blocks contribute comparably. Per-view scaling divides a whole block by its total variance/first singular value; reserve per-feature unit-variance scaling for inside PLS and only after filtering low-variance features.

```r
drop_constant <- function(mat, min_sd=1e-8) mat[apply(mat, 1, sd, na.rm=TRUE) > min_sd, ]   # per-feature scaling blows up zero-variance features

scale_per_view <- function(mat) mat / sqrt(sum(apply(mat, 1, var, na.rm=TRUE)))             # whole-block scaling: equalizes contribution, preserves within-block feature ratios

blocks <- lapply(list(RNA=vst_rna, Protein=norm_prot, Methylation=m_values), drop_constant)
blocks <- lapply(blocks, scale_per_view)
sapply(blocks, function(x) sum(apply(x, 1, var)))    # each block now contributes comparably
```

## Cross-Omic Batch Strategy

**Goal:** Remove technical batch once, in one place, without deleting biology or double-correcting.

**Approach:** First cross-tabulate batch against the biological variable; if they are collinear, stop - the effect is unrecoverable. Otherwise correct PER OMIC (ComBat on transformed data, ComBat-seq on counts) when the integrator needs clean input, OR model batch as a covariate inside the downstream step - never both, and never on a stacked multi-omic matrix.

```r
library(sva)

with(as.data.frame(colData(mae)), table(Batch, Condition))   # confounding gate: any empty cell = collinear -> do NOT correct

mod <- model.matrix(~ Condition, data=as.data.frame(colData(mae)))   # protect biology
vst_rna_bc <- ComBat(dat=vst_rna, batch=colData(mae)$Batch, mod=mod, par.prior=TRUE)   # ONE omic at a time
```

Stacking RNA, protein, and methylation into one matrix and running a single ComBat (with omic-type as the batch or a covariate) is a failure mode, not a recipe: it treats cross-omic differences as noise and erases the very signal the integration is meant to find. MOFA's multi-group `group=` is also not batch correction - it asks whether the same factors operate within each group, and does not regress batch out.

## Missing-Value Triage

**Goal:** Match the imputation to the missingness mechanism, and never fabricate a whole assay.

**Approach:** Separate the three regimes - missing features (filter), missing values within a feature (impute by mechanism: MNAR below detection vs MAR sporadic), and whole missing samples (do not impute; use a missing-tolerant integrator). Proteomics/metabolomics missingness is largely MNAR (a peptide is absent because it is low), so MAR methods bias it upward.

```r
keep <- rowMeans(is.na(norm_prot)) < 0.30          # drop features missing in >30% of samples before imputing
prot_f <- norm_prot[keep, ]

library(imputeLCMD)
prot_mnar <- impute.QRILC(prot_f)[[1]]             # left-censored draw for below-detection (MNAR) gaps
```

For a mosaic cohort (a subject profiled for RNA and methylation but not proteomics), do not impute the missing proteomics profile - MOFA2 ignores missing entries in its likelihood and tolerates incomplete views natively, so route that case to mofa-integration rather than fabricating a block that would manufacture cross-omic correlation.

## Per-Method Failure Modes

### Mismatched scales fed to a shared-latent method
**Trigger:** stacking raw counts, betas, and z-scores into one factorization. **Mechanism:** the model assumes comparable, homoscedastic, roughly-Gaussian features. **Symptom:** the leading factors are dominated by the highest-dynamic-range block; small blocks never surface. **Fix:** per-omic transform then per-view scaling, in that order.

### Per-feature scaling inflating noise
**Trigger:** mixOmics `scale=TRUE` (its default) on unfiltered blocks. **Mechanism:** a near-constant feature's tiny SD is divided out, blowing it up to unit variance. **Symptom:** a factor built from technical jitter at the detection floor. **Fix:** filter near-constant features before any per-feature scaling; prefer per-view scaling for cross-block equalization.

### Correcting confounded batch
**Trigger:** ComBat with batch correlated with (or equal to) the condition. **Mechanism:** the batch and biology terms are collinear, so the correction has no way to separate them. **Symptom:** the condition effect vanishes after correction. **Fix:** cross-tabulate batch x condition first; if collinear, do not correct - the design cannot be rescued post hoc.

### Double-correcting batch
**Trigger:** ComBat per omic AND a batch term inside the integration/DE step. **Mechanism:** the second correction operates on residuals already partly stripped. **Symptom:** over-shrunk effects and inflated confidence. **Fix:** correct once, in one place; prefer modeling for inferential steps, scrubbing only when the tool cannot accept a covariate.

### MAR imputation on MNAR data
**Trigger:** kNN/missForest on below-detection proteomics gaps. **Mechanism:** MAR methods borrow toward the observed mean, but the missing values are low by definition. **Symptom:** low-abundance proteins biased upward; the low-abundance biology is erased. **Fix:** QRILC/MinProb for MNAR entries, kNN/missForest only for sporadic MAR gaps.

### Feature filter treated as housekeeping
**Trigger:** top-N HVG / variable-feature selection per omic applied without justification. **Mechanism:** an integration model can only place a feature on a latent factor if the feature is in the input, so the cutoff pre-commits which axes of variation can be discovered. **Symptom:** two analysts with different cutoffs fit different models and get different factors; a real low-variance signal silently cannot appear. **Fix:** treat the filter as a modeling choice - declare and justify it, filter the wider views harder to equalize dimensionality, and confirm the headline factors are stable across a sensible cutoff range.

### Whole-sample imputation to satisfy a complete-case method
**Trigger:** imputing an entire missing omic profile so DIABLO/CCA will run. **Mechanism:** the imputed block is reconstructed from the other omics. **Symptom:** "discovered" cross-omic correlation that the imputation manufactured. **Fix:** use MOFA2 (models the missingness) or restrict to overlapping samples; never fabricate a whole assay.

### Orientation transpose and one-to-many ID mapping
**Trigger:** exporting an MAE assay to Python, or joining omics on gene symbols. **Mechanism:** Bioconductor is samples-in-columns while AnnData/mixOmics are samples-in-rows; symbols are non-unique and IDs map one-to-many. **Symptom:** genes treated as samples (silent), or duplicated/lost rows after a merge. **Fix:** assert matrix shape after every cross-language hop; join on Ensembl/UniProt/RefMet, and decide the collapse rule for ambiguous mappings explicitly.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Per-view scaling (not per-feature) for cross-block equalization | Escofier 1994 *Comput Stat Data Anal* 18:121 (MFA) | equalizes block contribution without inflating individual noise features |
| Near-constant feature filter (sd > ~1e-8) before per-feature scaling | mixOmics docs | per-feature unit-variance scaling blows up zero-variance features |
| Feature missingness filter ~30-50% before imputing | Lazar 2016 *J Proteome Res* 15:1116 | a feature missing in most samples cannot be imputed reliably |
| Batch x condition cross-tab must have no empty cell before ComBat | Nygaard 2016 *Biostatistics* 17:29 | an empty cell means batch and biology are collinear and inseparable |
| Model batch as covariate (not scrub) for inferential steps | Nygaard 2016 *Biostatistics* 17:29 | scrub-then-test understates variance and exaggerates significance |
| M-values (not beta) for methylation modeling | Du 2010 *BMC Bioinformatics* 11:587 | beta is heteroscedastic; the logit is approximately homoscedastic |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| One omic dominates every shared factor | no per-view scaling / unequal feature counts | per-view scale; filter the wider view harder |
| A factor is built from low-signal features | per-feature scaling without filtering | drop near-constant features before scaling |
| Condition effect disappears after batch correction | batch confounded with condition | check the cross-tab; do not correct collinear designs |
| Low-abundance proteins look unexpectedly high | MAR imputation on MNAR gaps | QRILC/MinProb for below-detection values |
| Cross-omic correlation that does not replicate | whole-sample imputation | use MOFA missing-view handling; do not fabricate a block |
| Genes appear as samples after export | orientation flip across languages | transpose; assert shape after every hop |

## References

- Ramos M, Schiffer L, Re A, et al. 2017. Software for the integration of multiomics experiments in Bioconductor. *Cancer Res* 77:e39-e42.
- Du P, Zhang X, Huang C-C, et al. 2010. Comparison of Beta-value and M-value methods for quantifying methylation levels by microarray analysis. *BMC Bioinformatics* 11:587.
- Love MI, Huber W, Anders S. 2014. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol* 15:550.
- Gloor GB, Macklaim JM, Pawlowsky-Glahn V, Egozcue JJ. 2017. Microbiome datasets are compositional: and this is not optional. *Front Microbiol* 8:2224.
- Escofier B, Pages J. 1994. Multiple factor analysis (AFMULT package). *Comput Stat Data Anal* 18:121-140.
- Johnson WE, Li C, Rabinovic A. 2007. Adjusting batch effects in microarray expression data using empirical Bayes methods. *Biostatistics* 8:118-127.
- Zhang Y, Parmigiani G, Johnson WE. 2020. ComBat-seq: batch effect adjustment for RNA-seq count data. *NAR Genom Bioinform* 2:lqaa078.
- Ritchie ME, Phipson B, Wu D, et al. 2015. limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Res* 43:e47.
- Leek JT, Storey JD. 2007. Capturing heterogeneity in gene expression studies by surrogate variable analysis. *PLoS Genet* 3:e161.
- Risso D, Ngai J, Speed TP, Dudoit S. 2014. Normalization of RNA-seq data using factor analysis of control genes or samples. *Nat Biotechnol* 32:896-902.
- Lazar C, Gatto L, Ferro M, Bruley C, Burger T. 2016. Accounting for the multiple natures of missing values in label-free quantitative proteomics data sets to compare imputation strategies. *J Proteome Res* 15:1116-1125.
- Nygaard V, Rodland EA, Hovig E. 2016. Methods that remove batch effects while retaining group differences may lead to exaggerated confidence in downstream analyses. *Biostatistics* 17:29-39.

## Related Skills

- integration-design - The method-selection decision this harmonization feeds
- mofa-integration - Consumes harmonized blocks; models missing-view samples natively
- mixomics-analysis - Consumes harmonized blocks; needs complete cases and is per-feature scaled
- similarity-network - Consumes harmonized blocks for patient stratification
- differential-expression/batch-correction - Single-omic RNA-seq batch mechanics (ComBat/limma)
- methylation-analysis/array-preprocessing - Methylation beta/M-value normalization
- proteomics/proteomics-qc - Proteomics normalization and QC
- metabolomics/normalization-qc - Metabolomics normalization and scaling
- metagenomics/abundance-estimation - Compositional/CLR theory for compositional omics

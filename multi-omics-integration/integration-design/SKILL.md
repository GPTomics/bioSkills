---
name: bio-multi-omics-integration-design
description: Chooses a bulk multi-omics integration strategy before any tool runs by mapping the biological question (subtype discovery, shared axis of variation, predictive signature, pairwise correlation) to a method class, naming the sample correspondence (paired-vertical, horizontal, mosaic, diagonal), enforcing the n<<p discipline that makes a held-out cohort the endpoint instead of in-cohort cross-validation, and running the per-view variance-imbalance diagnostic. Covers the early/mixed/intermediate/late taxonomy, why vertical and horizontal integration are different problems, and why a shared factor dominated by one omic is not integration. Use when deciding which integration method fits a question, whether data is paired or mosaic, supervised or unsupervised, or how to validate an integrated result. For unsupervised factors see mofa-integration; for supervised signatures see mixomics-analysis; for stratification see similarity-network; for single-cell see single-cell/multimodal-integration.
tool_type: r
primary_tool: MultiAssayExperiment
---

## Version Compatibility

Reference examples tested with: MultiAssayExperiment 1.36+, SummarizedExperiment 1.40+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

The tool versions that matter most are MOFA2 and mixOmics, whose APIs have moved across releases; this skill routes to those tool skills rather than calling them, so the binding version here is MultiAssayExperiment (the container in which the paired-vs-mosaic decision is made).

# Multi-Omics Integration Design

**"How should I integrate these omics?"** -> Map the biological question and the sample correspondence to a method class BEFORE running a tool - because at tens of samples and 10^5-10^6 features a spurious cross-omic signal is the default outcome, not the surprise.
- R: assemble a `MultiAssayExperiment`, then choose MOFA2 (shared factors) / mixOmics (signature) / SNF (subtypes) by question

Scope: the integration decision itself - method selection, correspondence (paired/horizontal/mosaic/diagonal), supervised-vs-unsupervised mapping, the n<<p discipline, and the variance-imbalance diagnostic. Running the chosen tool -> mofa-integration, mixomics-analysis, similarity-network. Cross-omic preprocessing -> data-harmonization. Single-cell multimodal -> single-cell/multimodal-integration. Horizontal same-feature meta-analysis -> differential-expression/batch-correction.

## The Single Most Important Modern Insight -- Bulk Multi-Omics Is a Small-n, Huge-p Discovery Problem Where a Spurious Cross-Omic Signal Is the Default

A typical bulk cohort has n = 30-300 samples and 10^4-10^6 features per omic, so after stacking blocks n is smaller than p by three to four orders of magnitude. In that regime an integrated signature that has not been validated out-of-sample is overwhelmingly noise that fit the training samples. The deliverable is never "the integrated signature" - it is question-matched structure that survives three gates, each of which a common failure violates:

1. **The correspondence gate.** Vertical integration (different omics, SAME samples) and horizontal integration (same features, different cohorts) are different problems. The bulk joint-latent tools are indexed by sample; feed them unpaired same-feature data and they still run but emit factors that are pure batch. Name the correspondence first.
2. **The variance gate.** Total variance scales with feature count and feature scale, so a 850k-CpG block out-votes a 100-metabolite block and the shared factors become methylation PCs. Inspect the per-factor, per-view variance-explained table every time; if one view dominates every shared factor, the integration re-discovered the biggest omic.
3. **The validation gate.** With n = 40 and 5-fold CV each fold tests 8 samples, so in-cohort cross-validation is optimistically biased to the point of fiction. An integrated subtype or signature is not credible until it reproduces in an INDEPENDENT cohort. The held-out cohort is the finding.

Organize the analysis around defending these three gates, not around picking a favorite tool.

## The Integration Taxonomy -- Three Orthogonal Axes, Not One Label

A method is a point in a 3D space, not a single name. Stating where a method sits on each axis prevents the category's two deepest errors (horizontal/vertical confusion and concatenation at n<<p).

| Axis | Values | What it decides |
|------|--------|-----------------|
| Stage - WHEN blocks combine (Ritchie 2015, Picard 2021) | early (concatenate then model), mixed (transform each block then combine), intermediate (jointly model blocks into shared + specific factors), late (model each omic, combine results) | early is worst at n<<p and variance imbalance; intermediate joint-latent (MOFA/JIVE/iCluster) is the discovery sweet spot; late is robust to missing blocks but drops feature-level cross-talk |
| Correspondence - WHAT is tied together (Argelaguet 2021) | vertical (diff omics, same samples - THIS category), horizontal (same features, diff cohorts - meta-analysis), mosaic (partial overlap), diagonal (no shared axis - single-cell) | conflating horizontal and vertical is the deepest category error; only vertical and mosaic belong here |
| Supervision - WHETHER an outcome drives it | unsupervised (discover subtypes/factors), supervised (predict/discriminate a label) | unsupervised plus then-correlate-with-outcome is hypothesis-generating, NOT a validated predictor |

MOFA = intermediate, vertical, unsupervised. DIABLO = intermediate, vertical, supervised. SNF = mixed/transformation, vertical, unsupervised. ComBat-across-cohorts = horizontal, unsupervised harmonization (routes OUT to differential-expression/batch-correction).

## Tool Taxonomy

| Tool / class | Citation | Stage / supervision | When |
|--------------|----------|---------------------|------|
| MOFA2 | Argelaguet 2018 *Mol Syst Biol* 14:e8124; Argelaguet 2020 *Genome Biol* 21:111 | intermediate, unsupervised | shared vs view-specific factors; tolerant of missing omics-per-sample; the default factor model -> mofa-integration |
| mixOmics DIABLO (`block.splsda`) | Singh 2019 *Bioinformatics* 35:3055; Rohart 2017 *PLoS Comput Biol* 13:e1005752 | intermediate, supervised | sparse cross-omic signature that DISCRIMINATES known groups -> mixomics-analysis |
| mixOmics sPLS (`spls`) | Rohart 2017 *PLoS Comput Biol* 13:e1005752 | intermediate, unsupervised | covariance-maximizing feature pairs between TWO blocks -> mixomics-analysis |
| mixOmics MINT (`mint.splsda`) | Rohart 2017 *BMC Bioinformatics* 18:128 | horizontal | SAME omic across multiple STUDIES (study as a known effect) - not cross-omic |
| SNF (SNFtool) | Wang 2014 *Nat Methods* 11:333 | mixed/transformation, unsupervised | patient stratification; feature count buys no votes; robust as complexity grows -> similarity-network |
| iCluster / iClusterPlus / moCluster | Shen 2009 *Bioinformatics* 25:2906; Meng 2016 *J Proteome Res* 15:755 | intermediate, unsupervised | ONE joint-latent clustering (vs reconciling K separate clusterings); subtype discovery |
| JIVE | Lock 2013 *Ann Appl Stat* 7:523 | intermediate, unsupervised | explicit joint + individual + noise decomposition (how much signal is cross-omic) |
| MFA / mixKernel | Mariette 2018 *Bioinformatics* 34:1009 | mixed | block weighting / kernel fusion to stop one omic dominating |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Different omics on the SAME samples, no phenotype, find shared axes | MOFA2 | unsupervised factor model; variance decomposition; native missing-block handling -> mofa-integration |
| Different omics on the same samples, want patient SUBTYPES | SNF + spectral clustering (or iCluster) | transformation-stage; robust to high p and a noisy omic -> similarity-network |
| Have a class label, want a cross-omic signature that discriminates it | mixOmics DIABLO + held-out cohort | supervised sparse multi-block PLS-DA -> mixomics-analysis |
| Just two omics, want correlated feature pairs | mixOmics sPLS | sparse PLS for a block pair -> mixomics-analysis |
| Quantify how much variation is joint vs omic-specific | JIVE (or the MOFA variance table) | explicit joint/individual split |
| SAME omic across multiple studies/cohorts | -> differential-expression/batch-correction or mixOmics MINT | horizontal integration / meta-analysis, NOT cross-omic |
| Mosaic cohort (some samples missing an omic) | MOFA2 (models the missingness) | intersecting to complete cases wastes scarce n -> data-harmonization |
| Single-cell CITE-seq / 10x Multiome / unpaired diagonal | -> single-cell/multimodal-integration | per-cell generative models; n is large; different paradigm |
| Per-omic DE then overlap the hit lists | -> differential-expression, methylation-analysis, proteomics | that is late integration by intersection, not joint modeling |
| Validate a discovered subtype against outcome | -> clinical-biostatistics/survival-analysis | survival / KM / Cox lives there |

Default when uncertain: assemble a `MultiAssayExperiment`, confirm vertical paired (or mosaic) correspondence, run MOFA2 for an unsupervised map and read its per-view variance-explained table, then escalate to a supervised (DIABLO) or stratification (SNF) tool only if the question demands it.

## Name the Correspondence First

**Goal:** Decide whether the data is a job for this category at all, and whether to model the missingness or intersect to complete cases.

**Approach:** Assemble the blocks into a `MultiAssayExperiment` (it coordinates assays, a sample map, and colData), then read off whether samples are fully paired, mosaic, or actually horizontal. Only vertical-paired and mosaic belong here.

```r
library(MultiAssayExperiment)

mae <- MultiAssayExperiment(experiments=ExperimentList(rna=rna_mat, prot=prot_mat, methyl=methyl_mat),
                            colData=clinical)
upsetSamples(mae)                 # visualize which samples have which omics (mosaic structure)
table(complete.cases(mae))        # how many samples have EVERY omic
paired <- intersectColumns(mae)   # complete-case fallback - counts the n it would cost
```

If `complete.cases` keeps most samples, complete-case methods (mixOmics, SNF) are fine. If a large fraction is mosaic, prefer MOFA2 (it models missing-view samples in its likelihood) over intersecting, because at n<<p discarding incomplete samples is expensive and imputing a whole block fabricates data (data-harmonization owns that decision).

## The Variance-Imbalance Diagnostic

**Goal:** Detect, before trusting any shared factor, whether one omic is set to dominate the integration purely because it has more features or larger scale.

**Approach:** After per-feature scaling, compare each block's total variance and feature count; a block contributing the overwhelming majority of stacked variance will hijack the shared latent space. The definitive check is the per-view variance-explained table that MOFA2 reports after fitting - if every factor loads on one view, equalize the blocks (MFA weighting, per-block keepX, or move to SNF) and refit.

```r
block_var <- sapply(assays_list, function(x) sum(apply(x, 1, var)))   # total variance per block
share     <- block_var / sum(block_var)
share                                                                 # any block >> others = imbalance risk
```

A block holding most of the stacked variance is a red flag that concatenation-style integration will re-discover it. This is the single best honesty check in the category; never skip the post-fit per-view variance read-out.

## The n<<p Discipline

The held-out cohort is the endpoint, not in-cohort cross-validation. Three rules follow from n<<p:

- **Tune with repeated cross-validation, never a single run.** At n = 40 a single CV estimate is mostly noise; mixOmics `perf`/`tune.*` take `nrepeat` (10-50) - use it. Generic CV/overfitting theory lives in machine-learning/model-validation.
- **Report out-of-sample performance, not the in-sample fit.** A supervised signature's training-set discrimination is guaranteed by construction; only an independent cohort makes it a biomarker.
- **An unsupervised factor that correlates with the outcome is a hypothesis.** MOFA found it without the label, which is a strength - but calling it predictive requires held-out validation, not the in-cohort correlation that found it.
- **Prefer regularized/sparse methods over early concatenation plus plain CCA.** Classical CCA divides out within-block variance and, when p > n, achieves correlation 1 trivially by overfitting; sparse PLS (covariance plus an L1 penalty) and sparse factor models stay identifiable, so the regularization is what makes the fit real, not a stylistic choice.

## Per-Method Failure Modes

### Horizontal data fed to a vertical method
**Trigger:** running MOFA/DIABLO/SNF on same-feature, multi-cohort data ("integrate my three RNA-seq studies"). **Mechanism:** the shared latent is indexed by sample and has nothing to align across feature-identical cohorts. **Symptom:** the tool runs and the top factors track cohort/run, not biology. **Fix:** recognize this as horizontal integration; use MINT, ComBat/sva, or differential-expression/batch-correction.

### Unvalidated in-cohort signature reported as a result
**Trigger:** reporting a DIABLO panel or a MOFA-factor-vs-outcome correlation from one cohort. **Mechanism:** at n<<p thousands of cross-omic feature pairs clear any threshold under the null; in-cohort CV is optimistically biased. **Symptom:** a beautiful signature that fails to replicate. **Fix:** hold out an independent cohort; frame an unvalidated finding as hypothesis-generating, never as a biomarker.

### Variance imbalance mistaken for integration
**Trigger:** concatenating blocks of very different feature counts/scales without equalization. **Mechanism:** the high-feature/high-variance omic casts the most votes for the shared factors. **Symptom:** every shared factor loads almost entirely on one view. **Fix:** read the per-view variance-explained table; equalize via MFA weighting / per-block keepX / per-feature z-scoring, or use SNF where each omic is one n x n network.

### Question-method mismatch
**Trigger:** using a tool whose output does not answer the question (e.g. SNF clusters reported with "driver features"). **Mechanism:** SNF selects no features, MOFA is unsupervised, DIABLO needs a label. **Symptom:** claims the method cannot support (SNF drivers without a post-hoc per-omic test; MOFA factors called predictive). **Fix:** map question -> class first (decision tree); do post-hoc per-omic differential analysis to find SNF subtype drivers.

### Cross-omic batch masquerading as shared biology
**Trigger:** omics generated on different platforms/labs/dates, interpreted without a technical check. **Mechanism:** the samples that ran together in every assay form a shared technical axis. **Symptom:** the top shared factor tracks run date / plate / site better than phenotype. **Fix:** correlate top factors against technical covariates before interpreting; correct per omic or model batch as a covariate (data-harmonization), watching for over-correction.

### Mosaic cohort forced to complete cases
**Trigger:** `intersectColumns` on a mosaic cohort before integrating. **Mechanism:** complete-case intersection drops every sample missing any omic. **Symptom:** n halves and power collapses. **Fix:** prefer MOFA2's native missing-view handling; reserve intersection for when mosaicism is minor.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| n<<p by ~3-4 orders of magnitude is the default regime | Subramanian 2020 *Bioinform Biol Insights* 14 | tens of samples, 10^4-10^6 features; dictates regularized/sparse methods and held-out validation |
| Repeated CV `nrepeat` 10-50 for any tuning at small n | mixOmics docs; n<<p variance | a single CV run at n~40 is noise; repetition stabilizes the estimate |
| Per-view variance-explained dominance flag: one view >~80% of every shared factor | Argelaguet 2018 *Mol Syst Biol* 14:e8124 (per-view variance decomposition) | a factor dominated by one view is view-specific structure, not integration |
| Drop MOFA factors below ~1-2% variance explained in every view | Argelaguet 2018 *Mol Syst Biol* 14:e8124 | low-variance factors are noise/over-parameterization |
| Held-out independent cohort for any reported biomarker/subtype | Subramanian 2020 *Bioinform Biol Insights* 14 | in-cohort CV at n<<p is optimistically biased; replication is the endpoint |
| Cross-check the headline with a second method class | Cantini 2021 *Nat Commun* 12:124; Tini 2019 *Brief Bioinform* 20:1269; Pierre-Jean 2020 *Brief Bioinform* 21:2011 | no single best method and methods diverge; a real result should survive a second class (e.g. a factor model and a fusion clustering) |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Tool runs on multi-cohort data but factors are all batch | horizontal data in a vertical method | use MINT / ComBat; this is meta-analysis, not cross-omic integration |
| n drops sharply after assembling the object | complete-case intersection on a mosaic cohort | use MOFA2 missing-view handling; intersect only if mosaicism is minor |
| Shared factors explain mostly one omic | variance imbalance (feature count / scale) | equalize blocks (MFA / per-block keepX / z-score) or use SNF |
| Signature does not replicate in a new cohort | reported from in-cohort CV at n<<p | hold out an independent cohort before claiming a biomarker |
| Two methods give different subtypes | method-dependent result (benchmarks rank methods differently because each optimizes a different criterion - robustness vs clustering recovery vs feature selection) | report which method and why; cross-check; no universal best method |

## References

- Ritchie MD, Holzinger ER, Li R, Pendergrass SA, Kim D. 2015. Methods of integrating data to uncover genotype-phenotype interactions. *Nat Rev Genet* 16:85-97.
- Picard M, Scott-Boyer M-P, Bodein A, Perin O, Droit A. 2021. Integration strategies of multi-omics data for machine learning analysis. *Comput Struct Biotechnol J* 19:3735-3746.
- Argelaguet R, Cuomo ASE, Stegle O, Marioni JC. 2021. Computational principles and challenges in single-cell data integration. *Nat Biotechnol* 39:1202-1215.
- Argelaguet R, Velten B, Arnol D, et al. 2018. Multi-Omics Factor Analysis - a framework for unsupervised integration of multi-omics data sets. *Mol Syst Biol* 14:e8124.
- Argelaguet R, Arnol D, Bredikhin D, et al. 2020. MOFA+: a statistical framework for comprehensive integration of multi-modal single-cell data. *Genome Biol* 21:111.
- Shen R, Olshen AB, Ladanyi M. 2009. Integrative clustering of multiple genomic data types using a joint latent variable model with application to breast and lung cancer subtype analysis. *Bioinformatics* 25:2906-2912.
- Lock EF, Hoadley KA, Marron JS, Nobel AB. 2013. Joint and individual variation explained (JIVE) for integrated analysis of multiple data types. *Ann Appl Stat* 7:523-542.
- Meng C, Helm D, Frejno M, Kuster B. 2016. moCluster: identifying joint patterns across multiple omics data sets. *J Proteome Res* 15:755-765.
- Mariette J, Villa-Vialaneix N. 2018. Unsupervised multiple kernel learning for heterogeneous data integration. *Bioinformatics* 34:1009-1016.
- Rohart F, Gautier B, Singh A, Le Cao K-A. 2017. mixOmics: an R package for 'omics feature selection and multiple data integration. *PLoS Comput Biol* 13:e1005752.
- Singh A, Shannon CP, Gautier B, et al. 2019. DIABLO: an integrative approach for identifying key molecular drivers from multi-omics assays. *Bioinformatics* 35:3055-3062.
- Wang B, Mezlini AM, Demir F, et al. 2014. Similarity network fusion for aggregating data types on a genomic scale. *Nat Methods* 11:333-337.
- Tini G, Marchetti L, Priami C, Scott-Boyer M-P. 2019. Multi-omics integration - a comparison of unsupervised clustering methodologies. *Brief Bioinform* 20:1269-1279.
- Cantini L, Zakeri P, Hernandez C, et al. 2021. Benchmarking joint multi-omics dimensionality reduction approaches for the study of cancer. *Nat Commun* 12:124.
- Pierre-Jean M, Deleuze J-F, Le Floch E, Mauger F. 2020. Clustering and variable selection evaluation of 13 unsupervised methods for multi-omics data integration. *Brief Bioinform* 21:2011-2030.
- Subramanian I, Verma S, Kumar S, Jere A, Anamika K. 2020. Multi-omics data integration, interpretation, and its application. *Bioinform Biol Insights* 14:1177932219899051.

## Related Skills

- mofa-integration - Unsupervised shared-factor discovery (the default tool once correspondence is vertical)
- mixomics-analysis - Supervised DIABLO signatures, sPLS pairs, and MINT multi-study integration
- similarity-network - Patient stratification via similarity network fusion
- data-harmonization - Per-block normalization, scaling, batch, and the mosaic missing-omic decision
- single-cell/multimodal-integration - Single-cell CITE-seq/Multiome integration (different paradigm)
- differential-expression/batch-correction - Horizontal same-feature meta-analysis and batch correction
- machine-learning/model-validation - Cross-validation and overfitting theory for supervised integration
- clinical-biostatistics/survival-analysis - Survival validation of discovered subtypes
- pathway-analysis/gsea - Enrichment of integrated factor or signature features
- workflows/multi-omics-pipeline - End-to-end multi-omics integration pipeline

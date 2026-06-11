---
name: bio-multi-omics-mofa-integration
description: Discovers shared and view-specific latent factors across bulk multi-omics blocks (RNA-seq, proteomics, methylation) on a common sample axis with MOFA2's unsupervised Bayesian group factor model, then attributes per-view variance explained and interprets signed factor weights. Covers why a factor is an unsupervised axis of variance and not a pathway, why a factor that correlates with batch is a batch factor, why the per-view variance-explained table is the primary read-out rather than p-values, why raw counts in a Gaussian view make factor 1 the library-size factor, and why MOFA2 handles missing omics-per-sample natively. Use when integrating two or more bulk omics to find joint axes of variation, choosing factor count, labeling factors against metadata, or running enrichment on factor weights. For supervised discriminant integration see mixomics-analysis; for the method decision see integration-design; for single-cell see single-cell/multimodal-integration; for enrichment see pathway-analysis/gsea.
tool_type: r
primary_tool: MOFA2
---

## Version Compatibility

Reference examples tested with: MOFA2 1.12+ (Bioconductor), mofapy2 0.7+, muon 0.1+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('MOFA2')` then `?function_name` to verify parameters
- Python: `pip show mofapy2 muon` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

MOFA2 is R/Bioconductor and trains through the Python mofapy2 backend via basilisk/reticulate; the trained model serializes to an `.hdf5` file whose schema is version-coupled. The R defaults differ from the muon Python defaults (spikeslab_weights, ard_factors, seed 42 vs 1) - confirm per ecosystem.

# MOFA2 Integration

**"Find shared variation across my omics layers"** -> Learn unsupervised latent factors that decompose variation into shared and view-specific axes - because the factors are found WITHOUT the phenotype, so a factor need not separate the groups of interest, and the variance-explained table is the read-out.
- R: `create_mofa(data_list)` -> `prepare_mofa()` -> `run_mofa()`
- Python: `mofapy2` for training, `muon` / `mofax` for downstream

Scope: unsupervised cross-block factor modeling of bulk omics, variance decomposition, and signed-weight interpretation. Supervised discriminant integration -> mixomics-analysis. The method-selection decision -> integration-design. Per-omic normalization / HVG / transform-to-Gaussian -> data-harmonization. Enrichment mechanics -> pathway-analysis/gsea. Single-cell multimodal MOFA+ -> single-cell/multimodal-integration.

## The Single Most Important Modern Insight -- A MOFA Factor Is an Unsupervised Axis of Variance, Not a Pathway, and a Factor That Correlates With Batch Is a Batch Factor

MOFA is PCA generalized to multiple omics: it returns factor scores (Z), per-view loadings (W), and a variance decomposition (R^2 per factor per view) that says which axis is shared across modalities and which is view-specific. That decomposition is the output. Three properties every misuse forgets:

1. **Factors are unsupervised and blind to the phenotype.** MOFA never saw the groups, so "factor 3 separates cases from controls" is genuinely strong evidence - but only if K factors were not scanned to find the one that splits the groups. A factor becomes biology only after its weights are annotated, it correlates with a known covariate, and it is shown NOT to be a technical-covariate factor.
2. **MOFA greedily captures the largest variance.** An unregressed batch effect or a depth gradient is often the largest variance, so it becomes a top factor. Always run `correlate_factors_with_covariates` against batch/depth/plate and exclude technical factors from biological interpretation. A factor that tracks run date is a batch factor, not confounded biology.
3. **Factors are unordered and ARD-pruned.** Factor 1 is not "the most important" the way PC1 is. The honest deliverable is "factor k explains X% of variance in views {A,B}, its top weights enrich for pathway P at sign s, it correlates with covariate C and not with batch, and it recurs across seeds" - anything less is an axis dressed up as a finding.

## Tool Taxonomy (Unsupervised Integrators)

| Tool | Citation | Output | When |
|------|----------|--------|------|
| MOFA / MOFA+ | Argelaguet 2018 *Mol Syst Biol* 14:e8124; Argelaguet 2020 *Genome Biol* 21:111 | continuous factors + R^2 per (factor,view) | interpretable variance-attributed axes; missing-data tolerant; the default factor model |
| MEFISTO (within MOFA2) | Velten 2022 *Nat Methods* 19:179 | smooth factors along a covariate | samples carry a spatial or temporal coordinate |
| iCluster / iClusterPlus / iClusterBayes | Shen 2009 *Bioinformatics* 25:2906; Mo 2018 *Biostatistics* 19:71 | a hard sample clustering | discrete SUBTYPES (a partition) rather than axes |
| JIVE | Lock 2013 *Ann Appl Stat* 7:523 | explicit joint + per-view-individual + residual subspaces | quantify and separate joint vs dataset-specific structure |
| MCIA | Meng 2016 *J Proteome Res* 15:755 (moCluster); omicade4 | co-inertia ordination across views | fast exploratory co-structure / visual pathway exploration |
| mixOmics DIABLO | Singh 2019 *Bioinformatics* 35:3055 | supervised discriminant signature | the OUTCOME must drive the projection -> mixomics-analysis |

## Decision Tree by Deliverable

| Deliverable | Recommended | Why |
|-------------|-------------|-----|
| Interpretable axes attributed across views, hypothesis generation | MOFA / MOFA+ | continuous factors + variance decomposition; native missing-data |
| Discrete patient subtypes (a partition) | iCluster / iClusterPlus | the output IS a hard clustering -> integration-design, similarity-network |
| Explicitly separate joint from view-specific structure | JIVE | decomposition is joint + individual + residual by construction |
| Samples have a spatial/temporal coordinate | MEFISTO (within MOFA2) | GP prior gives smooth factors along the covariate |
| The outcome/class must drive the projection | -> mixomics-analysis (DIABLO) | supervised; the phenotype is in the model |
| Some samples missing a whole omic | MOFA (handles it natively) | the likelihood ignores missing entries; do not impute a block |
| Single-cell CITE-seq / Multiome MOFA+ | -> single-cell/multimodal-integration | single-cell object plumbing and stochastic inference |
| Enrich the factor weights | -> pathway-analysis/gsea | GSEA/ORA mechanics; here only `run_enrichment` on weights |

## Prepare the Views

**Goal:** Get each omic into the orientation and distribution MOFA assumes, so the factors reflect biology rather than measurement scale.

**Approach:** Each view must be features-by-samples (the transpose of the usual samples-by-features matrix), per-omic normalized and variance-stabilized upstream, and HVG-filtered so feature counts are within an order of magnitude across views. The per-omic transform and HVG selection are owned by data-harmonization.

```r
library(MOFA2)

common <- Reduce(intersect, list(colnames(rna), colnames(prot), colnames(meth)))   # shared samples; mosaic samples may be kept, see below
data_list <- list(RNA=rna[, common], Protein=prot[, common], Methylation=meth[, common])   # each features x samples, already transformed
mofa <- create_mofa(data_list)
plot_data_overview(mofa)        # shows views, samples, and the missing-data pattern (grey = missing, tolerated)
```

MOFA tolerates missing samples in a view (it ignores missing entries in the likelihood), so a mosaic cohort can be passed directly rather than intersected to complete cases - this is a core reason to choose MOFA when data is incomplete.

## Create and Train

**Goal:** Configure a factor model whose count and likelihoods match the data and the sample size, then train by variational inference.

**Approach:** Over-specify the factor count and let ARD prune, transform counts to a Gaussian likelihood rather than using Poisson, set a seed for reproducibility, and write the model to a versioned `.hdf5`.

```r
data_opts  <- get_default_data_options(mofa)       # scale_views=FALSE, center_groups=TRUE
model_opts <- get_default_model_options(mofa)      # num_factors=10, likelihoods='gaussian'
train_opts <- get_default_training_options(mofa)   # convergence_mode='fast', drop_factor_threshold=-1, stochastic=FALSE, seed=42

model_opts$num_factors <- 15                       # over-specify; ARD prunes inactive factors
model_opts$likelihoods <- c(RNA='gaussian', Protein='gaussian', Methylation='gaussian')   # transform counts upstream, prefer gaussian
train_opts$drop_factor_threshold <- 0.01           # drop factors explaining <1% variance in ALL views
data_opts$scale_views <- TRUE                      # equalize per-view variance if feature counts cannot be balanced by filtering

mofa <- prepare_mofa(mofa, data_options=data_opts, model_options=model_opts, training_options=train_opts)
mofa <- run_mofa(mofa, outfile=file.path(tempdir(), 'model.hdf5'), use_basilisk=TRUE)
```

## Read the Variance Decomposition (the Output)

**Goal:** Identify which factors are shared across views and which are view-specific before interpreting any of them.

**Approach:** The R^2 per factor per view is the central result: a factor active in two or more views is a shared axis, a factor active in one view is view-specific. Inspect this table first; factors with near-zero R^2 everywhere are noise.

```r
var_exp <- get_variance_explained(mofa)            # $r2_total, $r2_per_factor  <- the central output
plot_variance_explained(mofa)                      # heatmap: factors x views
plot_variance_explained(mofa, plot_total=TRUE)     # total variance explained per view
```

## Label the Factors (and Exclude Technical Ones)

**Goal:** Earn a biological label for a factor instead of asserting one from its existence.

**Approach:** Attach metadata after fitting (it is never used to train), correlate each factor with both biological and technical covariates, and exclude any factor that tracks batch/depth/plate from biological interpretation. Then run enrichment on the signed weights, treating the two poles of the axis separately.

```r
md <- metadata[unlist(samples_names(mofa)), ]
md$sample <- rownames(md)                       # samples_metadata<- requires a literal 'sample' column
samples_metadata(mofa) <- md
correlate_factors_with_covariates(mofa, covariates=c('condition', 'batch', 'depth'))   # a factor that correlates with batch IS a batch factor

# enrichment per sign - the two poles are biological opposites along one axis
up   <- run_enrichment(mofa, view='RNA', feature.sets=msig_binary_matrix, factors=1:5, sign='positive')
down <- run_enrichment(mofa, view='RNA', feature.sets=msig_binary_matrix, factors=1:5, sign='negative')
```

Multi-group MOFA (`group` in `create_mofa`) partitions samples so factor activity can differ across groups while weights stay shared - it asks "do the same axes operate within each group?" It is NOT batch correction and putting the phenotype in as a group does not make MOFA supervised. For samples with a spatial or temporal coordinate, MEFISTO (a GP-prior mode of MOFA2, via `mefisto_options` + `set_covariates`) learns factors that vary smoothly along that covariate.

## Per-Method Failure Modes

### Factor narrated as a mechanism
**Trigger:** "factor 1 represents immune activation" from the factor's existence. **Mechanism:** a factor is a direction of covariation the ARD prior kept; it has no intrinsic meaning. **Symptom:** a biological story with no enrichment, no covariate correlation, no replication. **Fix:** report the R^2 footprint, the signed-weight enrichment, and the known-covariate correlation; call it hypothesis-generating until validated.

### Batch factor mistaken for biology
**Trigger:** interpreting a top factor without a technical-covariate check. **Mechanism:** MOFA captures the largest variance, and unregressed batch is often largest. **Symptom:** the top factor tracks run date / plate better than phenotype. **Fix:** regress known batch out upstream (before HVG selection); always `correlate_factors_with_covariates` and exclude technical factors.

### Confirmation-bias factor scanning
**Trigger:** reporting the one of K factors that splits the groups. **Mechanism:** with enough factors one will split any grouping by chance. **Symptom:** a factor-phenotype association that does not replicate. **Fix:** pre-specify the test or correct for K; validate the chosen factor on a held-out cohort.

### Raw counts in a Gaussian view
**Trigger:** feeding un-transformed RNA counts to a gaussian likelihood. **Mechanism:** counts are heavy-tailed and mean-variance-coupled, so high-count genes carry the most raw variance. **Symptom:** factor 1 tracks library size / housekeeping genes. **Fix:** normalize and variance-stabilize per view upstream (data-harmonization), then use gaussian; reserve bernoulli for genuine binaries.

### Big modality eats the factors
**Trigger:** a 20k-gene view beside a 50-feature view, unequalized. **Mechanism:** bigger modalities are overrepresented in the factors. **Symptom:** every factor describes mostly the large view. **Fix:** HVG-filter to comparable feature counts and/or set `scale_views=TRUE` (which changes the R^2 interpretation to within-view relative variance).

### Too many factors at small n
**Trigger:** num_factors=20 with 30 samples. **Mechanism:** factor analysis needs sample size (the package floor is >15). **Symptom:** factors that fit noise and split the cohort by accident. **Fix:** request fewer factors, prune with `drop_factor_threshold`, confirm robustness across seeds, validate out-of-sample.

### Seed fragility unchecked
**Trigger:** building a story on one training run, especially with `stochastic=TRUE`. **Mechanism:** PCA init makes standard VI mostly reproducible, but local optima and stochastic inference still vary. **Symptom:** a headline factor that does not reappear on a retrain. **Fix:** set the seed AND retrain with a different seed/factor count; a robust axis recurs with factor-score correlation near 1.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| `num_factors` over-specified, ARD prunes | Argelaguet 2018 *Mol Syst Biol* 14:e8124 | a factor never allocated cannot be recovered; default cap is N-dependent (5 if N<=25, 15 if N<=1000) |
| `drop_factor_threshold` ~0.01 | MOFA2 docs | drop factors explaining <1% variance in ALL views; default -1 keeps all |
| Sample size floor > 15 | MOFA2 FAQ | factor analysis is only useful with adequate n; tens of samples overfit a generous factor count |
| `scale_views=TRUE` only when filtering cannot equalize | MOFA2 FAQ | bigger modalities are overrepresented; scaling equalizes per-view variance but changes R^2 reading |
| Transform counts to gaussian rather than poisson | MOFA2 FAQ | non-gaussian likelihoods are less-accurate approximations; transform if it can be defended |
| Robustness: factor recurs across seeds with |r|~1 | Argelaguet 2018 *Mol Syst Biol* 14:e8124 | a factor that does not reappear on retrain is fragile noise |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `create_mofa` orientation error or nonsense factors | views passed samples-by-features | transpose to features-by-samples |
| Factor 1 tracks sequencing depth | raw counts into a gaussian view | normalize + variance-stabilize per view first |
| Every factor describes one omic | variance imbalance | HVG-filter to comparable feature counts / `scale_views=TRUE` |
| Model will not converge / factors NaN | unscaled blocks or a constant feature | center/scale; drop zero-variance features upstream |
| A factor loads almost entirely on one sample | an outlier hijacking a factor | inspect and remove the outlier; refit |
| Headline factor vanishes on rerun | seed fragility / stochastic inference | set seed; confirm the factor recurs across retrains |

## References

- Argelaguet R, Velten B, Arnol D, et al. 2018. Multi-Omics Factor Analysis - a framework for unsupervised integration of multi-omics data sets. *Mol Syst Biol* 14:e8124.
- Argelaguet R, Arnol D, Bredikhin D, et al. 2020. MOFA+: a statistical framework for comprehensive integration of multi-modal single-cell data. *Genome Biol* 21:111.
- Velten B, Braunger JM, Argelaguet R, et al. 2022. Identifying temporal and spatial patterns of variation from multimodal data using MEFISTO. *Nat Methods* 19:179-186.
- Shen R, Olshen AB, Ladanyi M. 2009. Integrative clustering of multiple genomic data types using a joint latent variable model with application to breast and lung cancer subtype analysis. *Bioinformatics* 25:2906-2912.
- Lock EF, Hoadley KA, Marron JS, Nobel AB. 2013. Joint and individual variation explained (JIVE) for integrated analysis of multiple data types. *Ann Appl Stat* 7:523-542.
- Meng C, Helm D, Frejno M, Kuster B. 2016. moCluster: identifying joint patterns across multiple omics data sets. *J Proteome Res* 15:755-765.
- Singh A, Shannon CP, Gautier B, et al. 2019. DIABLO: an integrative approach for identifying key molecular drivers from multi-omics assays. *Bioinformatics* 35:3055-3062.
- Cantini L, Zakeri P, Hernandez C, et al. 2021. Benchmarking joint multi-omics dimensionality reduction approaches for the study of cancer. *Nat Commun* 12:124.

## Related Skills

- integration-design - The method-selection decision; MOFA is the default once correspondence is vertical
- mixomics-analysis - Supervised DIABLO/sPLS where the outcome drives the projection
- data-harmonization - Per-omic transform, HVG selection, and batch regression before MOFA
- similarity-network - Hard patient stratification alternative to soft factors
- single-cell/multimodal-integration - Single-cell MOFA+ (CITE-seq/Multiome) plumbing
- pathway-analysis/gsea - Enrichment of factor weights (mechanics)
- clinical-biostatistics/survival-analysis - Survival validation using factors as features
- workflows/multi-omics-pipeline - End-to-end multi-omics integration pipeline

---
name: bio-workflows-multi-omics-pipeline
description: Orchestrates VERTICAL bulk multi-omics integration (RNA + protein + methylation on the SAME samples) from harmonization to a validated result, routing to MOFA2 (shared factors), mixOmics/DIABLO (predictive signature), or SNF (patient subtypes). Use when confirming the correspondence is vertical (not horizontal same-features-different-cohorts), joining on a stable sample primary key rather than cbind on assumed row order, normalizing each omic in its OWN space and equalizing block variance BEFORE stacking (or the widest omic hijacks every shared factor), correcting batch ONCE in one place, and validating in a HELD-OUT cohort because in-cohort CV at n<<p is optimistically biased. Hands mechanism to the multi-omics-integration component skills; not a re-teach of any single step.
tool_type: r
primary_tool: MOFA2
workflow: true
depends_on:
  - multi-omics-integration/integration-design
  - multi-omics-integration/data-harmonization
  - multi-omics-integration/mofa-integration
  - multi-omics-integration/mixomics-analysis
  - multi-omics-integration/similarity-network
qc_checkpoints:
  - after_harmonization: "Sample overlap >80% across blocks; join on a stable primary key (MAE sampleMap / MuData obs), never cbind on row order"
  - after_per_block_norm: "Each omic normalized in its OWN space (VST/M-values/log2+MNAR/CLR); missingness per modality <20%"
  - after_scaling: "Per-VIEW variance equalized (scale_views) BEFORE stacking; no single view's feature count dominates"
  - after_integration: "Per-(factor,view) R2 balanced; drop factors <1-2% R2 in ALL views; no single view dominates every factor"
  - after_validation: "Held-out cohort (not in-cohort CV, biased at n<<p); batch correlated with every factor; batch corrected ONCE"
---

## Version Compatibility

Reference examples tested with: MOFA2 1.12+, mixOmics 6.26+, SNFtool 2.3+, clusterProfiler 4.10+, ggplot2 3.5+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Multi-omics Integration Pipeline

**"Integrate my multi-omics datasets"** -> Decide the strategy first, then orchestrate harmonization, the chosen integration method (MOFA2, mixOmics, or SNF), interpretation, and validation - because bulk multi-omics is small-n, huge-p, so an unvalidated integrated result is the default noise outcome.

Before any tool runs, settle the design with multi-omics-integration/integration-design: confirm the data is vertical (different omics on the SAME samples, not the same features across cohorts), map the question to a method (shared factors -> MOFA2, predictive signature -> DIABLO, patient subtypes -> SNF), and plan a held-out cohort because in-cohort cross-validation at small n is optimistically biased. Inspect the per-view variance-explained table to confirm no single omic dominates the shared structure.

## The governing principle

Bulk multi-omics is small-n, huge-p, so an UNVALIDATED integrated result is the DEFAULT noise outcome, not the exception. Every honesty decision is made at a seam before the integrator runs.

1. **The sample-matching key must be VERTICAL, and joined on a stable primary key.** Vertical = different features, SAME samples (what this category is FOR); horizontal (same features, different cohorts) is meta-analysis/batch, a different problem — conflating them is the deepest category error. The shared latent factor is indexed by SAMPLE, so establish the subject primary key ONCE and externalize it (MAE `sampleMap` / consistent MuData `obs_names`); NEVER `cbind`/`pd.concat` on assumed row order (silently mis-pairs subjects). Gene SYMBOLS are display labels, not join keys (MARCH1->MARCHF1, Excel corruption) — join on stable Ensembl IDs.
2. **Each omic is normalized in its OWN space FIRST, then block variance is equalized BEFORE stacking.** A shared-latent integrator decomposes total variance and assumes each feature is ~continuous/Gaussian; the non-negotiable per-omic transforms (RNA-seq -> VST/logCPM never raw counts; methylation -> M-values; prot/metab -> log2 + MNAR-aware; microbiome -> CLR) come first, then per-VIEW scaling (MOFA `scale_views=TRUE`) equalizes each block's CONTRIBUTION. Variance is additive across features, so a block with more features (850k CpGs vs 100 metabolites) casts more votes and hijacks the shared space regardless of biological importance. Both failures are silent.
3. **No factor/subtype/signature is credible until it reproduces in a HELD-OUT cohort.** In-cohort CV at n<<p is optimistically biased (tuning feature selection on the full data then reporting same-data CV inflates AUC by up to ~0.15); the gold standard is an external test set or fully nested CV. Spurious cross-omic correlation is the DEFAULT at p>>n. Batch is corrected ONCE, in ONE place (ComBat each omic AND modeling batch downstream removes it twice and over-shrinks real biology).

**The single best honesty check:** if every shared factor is dominated by one view, the pipeline did not integrate — it re-discovered the biggest omic.

## Made-once commitments

| Commitment | Consequence inherited downstream |
|------------|----------------------------------|
| Correspondence axis (VERTICAL, not horizontal) | Whether the integrator's math has anything to align on; conflating the two is the deepest category error |
| Sample primary key (externalized in MAE/MuData) | Correct subject-to-assay linkage; cbind-on-row-order silently mis-pairs subjects |
| Per-block normalization + per-view variance scaling | Whether the widest omic hijacks every factor; done before stacking, never after |
| Validation plan (held-out cohort) | Whether any factor/signature is credible; in-cohort CV at n<<p is optimistically biased |

## Pipeline Overview

```
RNA-seq Data ─────┐
                  │
Proteomics Data ──┼──> Data Harmonization ──> Integration ──> Factors/Components
                  │                                                    │
Metabolomics ─────┘                                                    ▼
                        ┌─────────────────────────────────────────────────────┐
                        │           multi-omics-pipeline                      │
                        ├─────────────────────────────────────────────────────┤
                        │  1. Data Preprocessing per Modality                 │
                        │  2. Sample Harmonization (matching samples)         │
                        │  3. Feature Selection/Filtering                     │
                        │  4. Integration (MOFA2 / mixOmics / SNF)            │
                        │  5. Factor/Component Interpretation                 │
                        │  6. Downstream Analysis                             │
                        └─────────────────────────────────────────────────────┘
                                                                       │
                                                                       ▼
                              Integrated Factors + Biomarker Signatures
```

## Complete MOFA2 Workflow

**Goal:** Discover unsupervised shared and view-specific factors across the harmonized omics, then interpret and validate them.

**Approach:** Harmonize to common samples, feature-select per view, train MOFA2, read the per-view variance decomposition first, label factors against biological and technical covariates, and run enrichment on the signed weights.

```r
library(MOFA2)
library(MOFAdata)
library(ggplot2)
library(tidyverse)

# === 1. LOAD AND HARMONIZE DATA ===
# RNA-seq data (samples x genes)
rna <- read.csv('rnaseq_normalized.csv', row.names = 1)
cat('RNA:', nrow(rna), 'samples,', ncol(rna), 'genes\n')

# Proteomics data (samples x proteins)
protein <- read.csv('proteomics_normalized.csv', row.names = 1)
cat('Protein:', nrow(protein), 'samples,', ncol(protein), 'proteins\n')

# Metabolomics data (samples x metabolites)
metab <- read.csv('metabolomics_normalized.csv', row.names = 1)
cat('Metabolites:', nrow(metab), 'samples,', ncol(metab), 'metabolites\n')

# Find common samples
common_samples <- Reduce(intersect, list(rownames(rna), rownames(protein), rownames(metab)))
cat('Common samples:', length(common_samples), '\n')

# Subset to common samples
rna <- rna[common_samples, ]
protein <- protein[common_samples, ]
metab <- metab[common_samples, ]

# === 2. FEATURE SELECTION ===
# Select most variable features per modality
select_variable <- function(data, n = 2000) {
    vars <- apply(data, 2, var, na.rm = TRUE)
    top_features <- names(sort(vars, decreasing = TRUE))[1:min(n, ncol(data))]
    data[, top_features]
}

rna_var <- select_variable(rna, n = 2000)
protein_var <- select_variable(protein, n = 1000)
metab_var <- select_variable(metab, n = 500)

# === 3. CREATE MOFA OBJECT ===
# Prepare data as list of matrices (features x samples)
data_list <- list(
    RNA = t(as.matrix(rna_var)),
    Protein = t(as.matrix(protein_var)),
    Metabolome = t(as.matrix(metab_var))
)

# Create MOFA object
mofa <- create_mofa(data_list)

# Add sample metadata (samples_metadata<- requires a literal 'sample' column)
sample_metadata <- read.csv('sample_metadata.csv')
rownames(sample_metadata) <- sample_metadata$sample_id
sample_metadata$sample <- sample_metadata$sample_id
samples_metadata(mofa) <- sample_metadata[common_samples, ]

# === 4. CONFIGURE AND TRAIN MODEL ===
# Data options
data_opts <- get_default_data_options(mofa)
data_opts$scale_views <- TRUE  # Scale each view

# Model options
model_opts <- get_default_model_options(mofa)
model_opts$num_factors <- 15  # Number of factors to learn

# Training options
train_opts <- get_default_training_options(mofa)
train_opts$maxiter <- 1000
train_opts$convergence_mode <- 'slow'
train_opts$seed <- 42

# Prepare and train
mofa <- prepare_mofa(mofa, data_options = data_opts,
                      model_options = model_opts,
                      training_options = train_opts)

cat('Training MOFA model...\n')
mofa <- run_mofa(mofa, outfile = 'mofa_model.hdf5', use_basilisk = TRUE)

# === 5. ANALYZE FACTORS ===
# Variance explained per factor per view
plot_variance_explained(mofa, max_r2 = 15)
ggsave('variance_explained.png', width = 10, height = 6)

# Factor values
factor_values <- get_factors(mofa)[[1]]

# Correlate factors with biological AND technical covariates; a factor that tracks batch is a batch factor
correlate_factors_with_covariates(mofa, covariates = c('condition', 'batch', 'depth'))

# Factor plots
plot_factor(mofa, factors = 1:4, color_by = 'condition', dot_size = 3)
ggsave('factor_scatter.png', width = 12, height = 10)

# === 6. INTERPRET FACTORS ===
# Get top weights per factor per view
for (f in 1:5) {
    cat('\nFactor', f, ':\n')
    weights <- get_weights(mofa, factors = f, as.data.frame = TRUE)

    for (view in unique(weights$view)) {
        view_weights <- weights[weights$view == view, ]
        view_weights <- view_weights[order(abs(view_weights$value), decreasing = TRUE), ]
        cat('  ', view, ':', paste(head(view_weights$feature, 5), collapse = ', '), '\n')
    }
}

# Heatmap of top features per factor
plot_top_weights(mofa, view = 'RNA', factors = 1:5, nfeatures = 10)
ggsave('top_weights_rna.png', width = 10, height = 8)

# === 7. ENRICHMENT ANALYSIS ===
library(clusterProfiler)
library(org.Hs.eg.db)

# Get RNA weights for factor 1
rna_weights <- get_weights(mofa, views = 'RNA', factors = 1)[[1]][, 1]
top_genes <- names(sort(abs(rna_weights), decreasing = TRUE))[1:200]

# GO enrichment -- use all RNA features as background (not the full genome)
all_rna_genes <- names(rna_weights)
ego <- enrichGO(gene = top_genes,
                universe = all_rna_genes,
                OrgDb = org.Hs.eg.db,
                keyType = 'SYMBOL',
                ont = 'BP',
                pvalueCutoff = 0.05)
ego <- simplify(ego, cutoff = 0.7, by = 'p.adjust')

dotplot(ego, showCategory = 15)
ggsave('factor1_enrichment.png', width = 8, height = 10)

# === 8. DOWNSTREAM: SURVIVAL ANALYSIS ===
library(survival)
library(survminer)

# Add factor values to metadata
surv_data <- data.frame(
    sample = rownames(factor_values),
    factor1 = factor_values[, 1],
    time = sample_metadata[rownames(factor_values), 'survival_time'],
    status = sample_metadata[rownames(factor_values), 'survival_status']
)

# Median split
surv_data$factor1_group <- ifelse(surv_data$factor1 > median(surv_data$factor1), 'High', 'Low')

# Kaplan-Meier
fit <- survfit(Surv(time, status) ~ factor1_group, data = surv_data)
ggsurvplot(fit, data = surv_data, pval = TRUE, risk.table = TRUE)
ggsave('survival_factor1.png', width = 8, height = 8)

# === 9. EXPORT RESULTS ===
# Factor values
write.csv(factor_values, 'mofa_factor_values.csv')

# Weights
all_weights <- get_weights(mofa, as.data.frame = TRUE)
write.csv(all_weights, 'mofa_weights.csv', row.names = FALSE)

cat('\nMOFA analysis complete!\n')
```

## mixOmics DIABLO Workflow

**Goal:** Build a supervised cross-omic signature that discriminates a known outcome, with an honest performance estimate.

**Approach:** Set the design matrix from the goal, tune the component count then keepX inside cross-validation folds with balanced error rate, fit, and report performance from data not used in tuning.

```r
library(mixOmics)

# === 1. PREPARE DATA ===
# Same preprocessing as above
X <- list(
    RNA = as.matrix(rna_var),
    Protein = as.matrix(protein_var),
    Metabolome = as.matrix(metab_var)
)

# Outcome variable
Y <- factor(sample_metadata[common_samples, 'condition'])

# === 2. DESIGN MATRIX (the central DIABLO decision) ===
# off-diagonal trades discrimination vs cross-block correlation: ~1 for a coherent network,
# <0.5 for prediction. 0.1 leans toward prediction and is tutorial convention, not a default.
design <- matrix(0.5, ncol = length(X), nrow = length(X),
                 dimnames = list(names(X), names(X)))
diag(design) <- 0

# === 3. TUNE MODEL ===
# Tune number of components
perf.diablo <- perf(block.splsda(X, Y, ncomp = 5, design = design),
                    validation = 'Mfold', folds = 5, nrepeat = 10)

ncomp <- perf.diablo$choice.ncomp$WeightedVote['Overall.BER', 'max.dist']
cat('Optimal components:', ncomp, '\n')

# Tune number of variables per component
test.keepX <- list(
    RNA = c(10, 25, 50, 100),
    Protein = c(5, 10, 25, 50),
    Metabolome = c(5, 10, 25)
)

tune.diablo <- tune.block.splsda(X, Y, ncomp = ncomp, test.keepX = test.keepX,
                                  design = design, validation = 'Mfold', folds = 5)

optimal.keepX <- tune.diablo$choice.keepX

# === 4. FINAL MODEL ===
diablo.model <- block.splsda(X, Y, ncomp = ncomp,
                              keepX = optimal.keepX, design = design)

# === 5. VISUALIZATION ===
# Sample plot
plotIndiv(diablo.model, ind.names = FALSE, legend = TRUE, title = 'DIABLO Sample Plot')

# Variable plot
plotVar(diablo.model, var.names = FALSE, style = 'graphics', legend = TRUE)

# Circos plot
circosPlot(diablo.model, cutoff = 0.7, line = TRUE,
           color.blocks = c('darkorchid', 'brown1', 'lightgreen'))

# Heatmap
cimDiablo(diablo.model, color.blocks = c('darkorchid', 'brown1', 'lightgreen'),
          margins = c(10, 5))

# === 6. PERFORMANCE (report from data not used to tune; an external test set is the honest estimate) ===
perf.final <- perf(diablo.model, validation = 'Mfold', folds = 5, nrepeat = 10)
perf.final$WeightedVote.error.rate   # matrix: classes + Overall.BER by component

# ROC curves
auc.diablo <- auroc(diablo.model, roc.block = 'RNA', roc.comp = 1)
```

## Similarity Network Fusion (SNF)

**Goal:** Stratify patients into candidate subtypes from the fused multi-omic similarity network.

**Approach:** Standardize each omic, build local-scaled affinity networks, fuse by cross-diffusion, estimate a plausible cluster number, and defend it with a fused-versus-single-omic concordance check before claiming subtypes.

```r
library(SNFtool)

# === 1. CREATE SIMILARITY MATRICES ===
K <- 20       # neighbors for the local kernel bandwidth (10-30)
sigma <- 0.5  # affinityMatrix width (the arg is sigma, not alpha); 0.3-0.8

# standardize per feature, then squared-Euclidean -> root -> local-scaled kernel
views <- lapply(list(rna_var, protein_var, metab_var), function(x) standardNormalization(as.matrix(x)))
affinities <- lapply(views, function(x) affinityMatrix(dist2(x, x)^(1/2), K, sigma))   # dist2 returns SQUARED distance

# === 2. FUSE NETWORKS ===
W <- SNF(affinities, K, t = 20)

# === 3. CLUSTER ON FUSED NETWORK (defend the count, do not assume it) ===
estimateNumberOfClustersGivenGraph(W, NUMC = 2:8)   # four eigengap/rotation estimates - plausibility, not truth
clusters <- spectralClustering(W, K = 3)            # here K is the CLUSTER COUNT
concordanceNetworkNMI(c(affinities, list(W)), 3)    # did fusion beat the best single omic? (Rappoport and Shamir 2018)

# === 4. VISUALIZATION ===
# Plot fused network
displayClustersWithHeatmap(W, clusters)
```

## QC Checkpoints

| Stage | Check | Action if Failed |
|-------|-------|------------------|
| Sample matching | >80% samples shared | Check sample IDs |
| Missing values | <20% per modality | Impute or remove |
| Feature variance | Features vary | Filter low variance |
| Model convergence | ELBO plateau | Increase iterations |
| Factor variance | drop factors below ~1-2% in all views | set drop_factor_threshold; keep fewer factors |
| Variance imbalance | no single view dominates every factor | per-view scaling or filter the wider view harder |
| Validation | held-out cohort, not in-cohort CV | replicate before claiming a biomarker/subtype |

## Workflow Variants

### With Missing Samples
```r
# MOFA2 handles missing views gracefully; create_mofa_from_df wants one row per (sample, feature, value)
to_long <- function(mat, view) {
    df <- as.data.frame(as.table(as.matrix(mat)))   # samples x features -> Var1=sample, Var2=feature, Freq=value (alignment preserved)
    data.frame(sample = as.character(df$Var1), feature = as.character(df$Var2), view = view, value = df$Freq)
}
data_long <- rbind(to_long(rna, 'RNA'), to_long(protein, 'Protein'))
mofa <- create_mofa_from_df(data_long)
```

### Single-cell Multi-omics
Single-cell multimodal data (CITE-seq, 10x Multiome) is a different paradigm - per-cell generative models with abundant observations rather than the bulk small-n regime. Route it to single-cell/multimodal-integration rather than applying this bulk pipeline.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Integrated result on mislabeled subjects (both axes have plausible lengths) | `cbind`/`pd.concat` on assumed sample order | Enforce sample linkage via a container (MAE `sampleMap` / MuData `intersect_obs`); assert shape after every cross-language hop |
| One view dominates every factor | Block scale/feature-count imbalance (850k CpGs vs 100 metabolites) | Per-VIEW scaling (`scale_views=TRUE`), NOT per-feature; filter larger views harder; check per-(factor,view) R2 |
| A near-constant noise feature blown up to variance 1 | Per-feature scaling (`scale=TRUE` in mixOmics) | Aggressive low-variance filtering BEFORE scaling; prefer per-view scaling |
| A "biological" factor that is really batch | Technical variance not regressed out | Correlate every factor with technical covariates (run date, plate, site); exclude a factor tracking batch more than biology |
| AUC inflated by up to ~0.15 | In-cohort CV with feature selection outside the folds | CV must WRAP selection (nested); gold standard is an external test set (`predict()`/`auroc()`) |
| Batch removed twice, real biology over-shrunk | ComBat each omic AND batch in the downstream model | Correct batch ONCE: pre-correct and do not re-enter, OR leave raw and model batch as a covariate |
| Transposed matrices make every gene a "sample" | Bioconductor (samples=cols) vs mixOmics/MOFA (samples=rows) mismatch | Transpose deliberately + assert shape on every cross-package hop |

## References

- Argelaguet R, Velten B, Arnol D, et al (2018) Multi-Omics Factor Analysis - a framework for unsupervised integration of multi-omics data sets. *Molecular Systems Biology* 14:e8124. DOI 10.15252/msb.20178124. (MOFA.)
- Rohart F, Gautier B, Singh A, Le Cao KA (2017) mixOmics: an R package for 'omics feature selection and multiple data integration. *PLoS Computational Biology* 13:e1005752. DOI 10.1371/journal.pcbi.1005752.
- Singh A, Shannon CP, Gautier B, et al (2019) DIABLO: an integrative approach for identifying key molecular drivers from multi-omics assays. *Bioinformatics* 35:3055-3062. DOI 10.1093/bioinformatics/bty1054.
- Wang B, Mezlini AM, Demir F, et al (2014) Similarity network fusion for aggregating data types on a genomic scale. *Nature Methods* 11:333-337. DOI 10.1038/nmeth.2810. (SNF.)
- Rappoport N, Shamir R (2018) Multi-omic and multi-view clustering algorithms: review and cancer benchmark. *Nucleic Acids Research* 46:10546-10562. DOI 10.1093/nar/gky889. (integration is not automatically better than the best single omic.)
- Nygaard V, Rodland EA, Hovig E (2016) Methods that remove batch effects while retaining group differences may lead to exaggerated confidence in downstream analyses. *Biostatistics* 17:29-39. DOI 10.1093/biostatistics/kxv027.

## Related Skills

- multi-omics-integration/integration-design - Method selection, correspondence, and the n<<p discipline (decide first)
- multi-omics-integration/mofa-integration - MOFA2 unsupervised factor analysis
- multi-omics-integration/mixomics-analysis - mixOmics DIABLO/sPLS/MINT methods
- multi-omics-integration/similarity-network - SNF patient stratification
- multi-omics-integration/data-harmonization - Cross-omic preprocessing and scaling
- pathway-analysis/go-enrichment - Factor/signature interpretation
- differential-expression/batch-correction - Batch effects
- clinical-biostatistics/survival-analysis - Survival validation of factors and subtypes
- single-cell/multimodal-integration - Single-cell multimodal integration (different paradigm)

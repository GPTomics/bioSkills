---
name: bio-single-cell-differential-abundance
description: Test whether cell-type proportions or composition changed between conditions in single-cell data using Milo (miloR), scCODA, sccomp, and propeller. Use when comparing cell-type proportions / composition between conditions, asking which populations expanded or contracted with treatment or disease, running neighborhood-level (cluster-free) abundance testing, or guarding against compositional shifts that masquerade as differential expression.
tool_type: mixed
primary_tool: Milo
---

## Version Compatibility

Reference examples tested with: miloR 2.0+, scCODA 0.1.9+, sccomp 1.8+, speckle 1.0+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Differential Abundance Testing

**"Did cell-type proportions change between conditions?"** -> Test whether populations expanded or contracted between groups, accounting for the fact that proportions are not independent.
- R (cluster-free): `miloR` - build kNN graph, define neighborhoods, `testNhoods()` with a GLM and SpatialFDR
- Python/R (cluster-based): `scCODA` (Bayesian Dirichlet-multinomial), `sccomp` (Bayesian, outlier-robust), `propeller` (speckle, arcsin-sqrt + limma)

## Governing principle

Composition data live on a SIMPLEX: proportions sum to 1, so they are NOT independent - when one population expands, every other proportion is mechanically forced down even if its absolute count never changed. Running a per-cluster t-test (or Wilcoxon) on proportions across samples is therefore invalid: it ignores the negative correlation the constraint imposes, treats each cell type as a free measurement, and produces correlated false positives (one true expansion drags down the rest, which then test as spurious "depletions"). Valid methods model the joint composition: either a Dirichlet-multinomial / log-ratio model with a reference (scCODA, sccomp) or a variance-stabilizing transform plus a linear model (propeller), or they sidestep hard clusters entirely by testing abundance on the kNN graph (Milo).
Replicates are samples, not cells. The unit of replication for a composition claim is the biological sample/donor; thousands of cells from one donor are one draw. Differential abundance needs biological replicates per condition (Milo, scCODA, sccomp, propeller all model sample-level counts), and few replicates (n<3-4/group) leave abundance shifts underpowered and unstable - more donors help, more cells per donor barely do. With n=1 per condition the donor is perfectly confounded with condition: the effect is unidentifiable, not merely underpowered, yet scCODA and sccomp will still emit confident `credible_effects()` that are pure donor idiosyncrasy - require >=2 (ideally 3-4) biological replicates per group before believing any abundance call.
Differential abundance and differential expression are different questions and confound each other. A pseudobulk or cluster-level "DE" signal between conditions can be pure composition: if a cluster mixes substates and treatment shifts their ratio, the aggregated profile changes although no gene changed expression in any cell - differential abundance masquerading as differential expression, invisible if only DE is run. Always pair a condition-DE analysis (single-cell/markers-annotation, differential-expression/deseq2-basics) with a differential-abundance test and interpret them jointly.

## Choosing a differential-abundance method

| Method | Model | Granularity | Use when | Fails when |
|--------|-------|-------------|----------|------------|
| Milo (miloR) | NB-GLM on kNN-neighborhood counts, SpatialFDR | Cluster-free neighborhoods | Continuous/transitional states; shifts that discrete clusters hide; want sub-cluster resolution | Very few cells/sample; results sensitive to k and `prop`; needs an integrated embedding |
| scCODA | Bayesian Dirichlet-multinomial, log-linear, reference cell type | Discrete clusters | Cluster-level testing with the simplex bias handled; want credible effects / FDR | Reference cell type mis-chosen; very few samples; HMC tuning |
| sccomp | Bayesian beta-binomial mixed model, outlier-robust | Discrete clusters | Outliers/over-dispersion present; want joint mean + variability, random effects | Small data with weak priors; longer runtime |
| propeller (speckle) | arcsin-sqrt or logit transform + limma moderated test | Discrete clusters | Fast frequentist test, several samples/group, Seurat/SCE input | Very small sample counts; ignores some compositional coupling vs Bayesian models |
| Simple proportion t-test / chi-square | Per-cluster test on proportions | Discrete clusters | Never recommended as the primary test | Always - ignores the simplex; correlated false positives |

scCODA and sccomp are cluster-based and Bayesian and report credible/FDR-controlled effects; Milo is cluster-free and catches shifts within a cell type that clustering averages away; propeller is the fast frequentist option. Run a cluster-based method and Milo when feasible and reconcile. When methods compete, verify current best practice against installed docs.

## The reference-cell-type choice in scCODA

Compositional analysis is always relative to something. scCODA fixes one cell type as the reference assumed unchanged by the covariates, and reports every other type's change relative to it; the verdict can flip with a different reference. Choose a cell type that is biologically stable and abundant across all samples, or use `reference_cell_type='automatic'` (scCODA picks a type with low dispersion present in all samples). A reference that actually changes will bias all other calls. sccomp avoids a hard reference by modeling all groups jointly; Milo avoids it via the graph.

## Adjusting for nuisance covariates and confounded designs

**Goal:** Adjust the abundance model for technical or biological nuisances (sequencing batch, timing, sex, age) and recognize when adjustment cannot help.

**Approach:** Add the nuisance as an extra additive term in the model formula with the condition of interest last; the test then reports the condition effect holding the nuisance constant. The nuisance column must vary within each condition - if a batch is perfectly confounded with condition (e.g. all controls sequenced in batch 1, all treated in batch 2), the term is unidentifiable and the test is invalid; the fix is experimental (multiplex conditions across batches), not statistical.

```r
# Milo: batch added before condition; batch column lives in design.df
design <- distinct(as.data.frame(colData(milo))[, c('sample', 'batch', 'condition')])
rownames(design) <- design$sample
da <- testNhoods(milo, design = ~ batch + condition, design.df = design, reduced.dim = 'PCA')
```

```python
# scCODA: additive patsy formula; covariate columns must be in the count table
data = dat.from_pandas(counts, covariate_columns=['sample', 'batch', 'condition'])
model = mod.CompositionalAnalysis(data, formula='batch + condition', reference_cell_type='automatic')
```

```r
# sccomp: nuisance added to formula_composition (and optionally formula_variability)
res <- sccomp_estimate(counts_tbl, formula_composition = ~ batch + condition, .sample = sample, .cell_group = cell_type, .count = count, cores = 1)
```

Diagnose confounding before modeling: cross-tabulate batch x condition; if a batch maps to a single condition, no covariate term recovers the effect. Build Milo's kNN graph on a batch-corrected embedding, but keep batch in the GLM design as well, since integration and design adjustment address different residual structure.

## Milo - cluster-free neighborhood abundance (R)

**Goal:** Test differential abundance on kNN neighborhoods so shifts within and between cell types are both visible.

**Approach:** Build the Milo object from an integrated reduced dimension, sample representative neighborhoods, count cells per sample per neighborhood, then fit a GLM with `testNhoods` and control the graph-aware SpatialFDR; annotate neighborhoods back to cell types for interpretation.

```r
library(miloR)
library(SingleCellExperiment)

milo <- Milo(sce)
milo <- buildGraph(milo, k = 30, d = 30, reduced.dim = 'PCA')
milo <- makeNhoods(milo, prop = 0.1, k = 30, d = 30, refined = TRUE, reduced_dims = 'PCA')
milo <- countCells(milo, meta.data = as.data.frame(colData(milo)), samples = 'sample')

design <- data.frame(colData(milo))[, c('sample', 'condition')]
design <- distinct(design)
rownames(design) <- design$sample
milo <- calcNhoodDistance(milo, d = 30, reduced.dim = 'PCA')

da <- testNhoods(milo, design = ~ condition, design.df = design, reduced.dim = 'PCA')
da <- annotateNhoods(milo, da, coldata_col = 'cell_type')
table(da$SpatialFDR < 0.1, da$cell_type)
```

`k` and `prop` trade resolution against power: larger neighborhoods are better powered but blur fine shifts. SpatialFDR (not raw p) corrects for overlapping neighborhoods - report it. A neighborhood with a mixed `cell_type` fraction is a genuinely transitional region, not a labeling error.

## scCODA - Bayesian cluster-level composition (Python)

**Goal:** Test cluster proportion changes while handling the simplex's negative-correlation bias.

**Approach:** Build a per-sample cell-type count table with covariates, fit the Dirichlet-multinomial model against a reference cell type, sample the posterior, then read credible effects at a chosen FDR.

```python
import pandas as pd
from sccoda.util import cell_composition_data as dat
from sccoda.util import comp_ana as mod

counts = pd.crosstab(adata.obs['sample'], adata.obs['cell_type']).reset_index()
meta = adata.obs[['sample', 'condition']].drop_duplicates()
counts = counts.merge(meta, on='sample')

data = dat.from_pandas(counts, covariate_columns=['sample', 'condition'])
model = mod.CompositionalAnalysis(data, formula='condition', reference_cell_type='automatic')
result = model.sample_hmc()
result.set_fdr(est_fdr=0.1)
result.summary()
print(result.credible_effects())
```

`set_fdr(est_fdr=0.1)` chooses the spike-and-slab threshold for the desired expected FDR; credible effects are the populations whose change is supported relative to the reference.

## sccomp - outlier-robust Bayesian composition (R)

**Goal:** Test composition (and variability) jointly, robust to outlier samples.

**Approach:** Estimate the beta-binomial model from a count table or cell-level data with `sccomp_estimate`, optionally remove outliers, then test contrasts with `sccomp_test`, which returns a Bayesian FDR (`c_FDR`).

```r
library(sccomp)

res <- counts_tbl |>
    sccomp_estimate(formula_composition = ~ condition, .sample = sample, .cell_group = cell_type, .count = count, cores = 1) |>
    sccomp_remove_outliers(cores = 1) |>
    sccomp_test()
res[res$c_FDR < 0.05, c('cell_type', 'c_effect', 'c_FDR')]
```

`sccomp_test` reports `c_effect` (composition log-fold change) and `c_FDR`; modeling variability separately catches groups that differ in dispersion, not just mean proportion.

## propeller - fast frequentist proportions (R)

**Goal:** Quickly test cell-type proportion differences across groups.

**Approach:** Compute per-sample proportions, apply an arcsin-sqrt (or logit) variance-stabilizing transform, and run a limma moderated test per cell type.

```r
library(speckle)

out <- propeller(clusters = seurat_obj$cell_type, sample = seurat_obj$sample, group = seurat_obj$condition)
out[out$FDR < 0.05, ]
```

propeller is the fast default for several samples per group; for outliers, over-dispersion, or random effects, prefer sccomp or scCODA.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Many cell types flagged as changed, all anti-correlated | Per-cluster proportion t-tests ignore the simplex | Use scCODA/sccomp/propeller/Milo, which model the joint composition |
| scCODA verdict flips between runs | Reference cell type mis-chosen or actually changing | Pick a stable abundant reference, or `reference_cell_type='automatic'` |
| No significant abundance change despite an obvious shift | Too few biological replicates; underpowered | Add donors (not cells); report effect sizes / credible intervals |
| Milo neighborhoods look noisy / unstable | k or `prop` too small, or embedding not integrated | Increase k/prop; build the graph on a batch-corrected reduced dim |
| "DE genes" between conditions but expression unchanged per cell | Compositional shift masquerading as DE | Run a differential-abundance test alongside the DE analysis |
| propeller p-values too liberal with few samples | Frequentist test under-powered/over-confident at small n | Use a Bayesian model (sccomp/scCODA) and report uncertainty |
| Abundance significant only in one direction across all types | Reporting raw proportions without the constraint | Interpret relative to a reference and report which population actually drives the shift |

## Related Skills

- clustering - Define the clusters whose abundance is tested (cluster-based methods)
- cell-annotation - Annotate cell types before testing their proportions
- markers-annotation - Pair condition DE with abundance testing to separate the confound
- batch-integration - Build the integrated embedding Milo's kNN graph relies on
- differential-expression/deseq2-basics - Pseudobulk condition DE that abundance testing complements
- pathway-analysis/go-enrichment - Characterize the populations that expanded or contracted

## References

- Dann et al. 2022, Nat Biotechnol 40:245-253 - Milo; differential abundance on kNN-graph neighborhoods with SpatialFDR.
- Buttner et al. 2021, Nat Commun 12:6876 - scCODA; Bayesian Dirichlet-multinomial compositional analysis with a reference cell type.
- Mangiola et al. 2023, PNAS 120(33):e2203828120 - sccomp; outlier-robust Bayesian differential composition and variability.
- Phipson et al. 2022, Bioinformatics 38(20):4720 - propeller; arcsin-sqrt transform plus limma for cell-type proportion testing.
- Squair et al. 2021, Nat Commun 12:5692 - sample, not cell, is the unit of replication for cross-condition single-cell claims.

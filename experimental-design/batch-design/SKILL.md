---
name: bio-experimental-design-batch-design
description: Designs genomics experiments so technical nuisance variation (batch, lane, plate, flow cell, operator, reagent lot, processing day) is balanced against the biological variable of interest and therefore estimable rather than confounded, using constrained sample-to-batch assignment (designit, OSAT), the confounder/mediator/collider distinction, and the principle that no post-hoc correction recovers a fully confounded design. Covers detecting hidden batches with surrogate variable analysis, a decision table for downstream correction (ComBat-seq, RUVSeq, SVA) whose execution is deferred to differential-expression/batch-correction, and reproducibility metadata. Use when assigning samples to sequencing batches/lanes/plates, avoiding batch-condition confounding, deciding whether a design is salvageable by correction, choosing a correction method, or estimating the number of hidden batches. For the experimental unit, randomization, and blocking concepts see experimental-design/randomization-blocking.
tool_type: r
primary_tool: designit
---

## Version Compatibility

Reference examples tested with: designit 0.5+, OSAT 1.50+ (Bioconductor), sva 3.50+, RUVSeq 1.36+, limma 3.58+, edgeR 4.0+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws an error, introspect the installed package and adapt to the actual API. Notes: OSAT uses `optimal.shuffle()` on a setup object (there is no bare `osat()` function); designit is R6 (`BatchContainer$new()`, `optimize_design()`, `*_score_generator()`) and its signatures drift between releases; `sva::ComBat_seq()` is for integer counts while `ComBat()` expects log-normalized values. Confirm each against the installed vignette before relying on it.

# Batch Design

**"Design my experiment so batch effects don't ruin it"** -> Assign samples to batches/lanes/plates so the biological variable is balanced against (orthogonal to) every technical nuisance factor, making batch estimable rather than confounded — because no post-hoc correction recovers a design where batch and condition are aliased.
- R: `designit::optimize_design()`, `OSAT::optimal.shuffle()` — constrained assignment at design time
- R: `sva::sva()`/`num.sv()` — detect hidden batches; `sva::ComBat_seq()`, `RUVSeq::RUVg()` — DOWNSTREAM correction (executed in differential-expression/batch-correction)

## The Single Most Important Modern Insight -- No Post-Hoc Method Recovers a Confounded Design

When a technical factor is perfectly aliased with the biological factor (all treated in batch 1, all controls in batch 2), batch and condition occupy the same column space and are mathematically non-identifiable; ComBat, SVA, and RUV cannot separate them, and "removing the batch effect" removes the biology with it. Worse, on *partially* confounded or merely **unbalanced** designs, mean-centering batches can *manufacture* false positives and inflate downstream confidence — Nygaard, Rødland & Hovig 2016 *Biostatistics* 17:29 showed a pipeline that returned >1000 spurious DE probes where the honest analysis (batch kept *in the model*) found 11. The operative rules: (1) balance the biological variable across batches **at design time** because correction is not a rescue (Leek 2010 *Nat Rev Genet* 11:733); (2) for inference, keep batch *in the model* so its degrees of freedom are charged honestly — reserve a batch-"cleaned" matrix for visualization and clustering only. The experimental unit, not the measurement, still defines replication (see experimental-design/randomization-blocking).

## Confounding vs Blocking vs Nuisance -- the Causal-Graph View

Batch is the genomics face of confounding, and not all metadata should be "adjusted for". A **confounder** is a common cause of treatment and outcome (adjust for it); a **mediator** lies on the causal path (adjusting removes signal); a **collider** is a common effect (adjusting *induces* spurious association). "Adjust for everything measured" is therefore wrong in general — conditioning on a collider opens a backdoor path. In a properly *randomized* design the biological variable has no confounders by construction, so batch is handled by **balanced assignment + a block/covariate term**, not by scrubbing every measured variable.

## Algorithmic Taxonomy -- Design Strategies for Technical Variation

| Strategy | What it does | When to use | Fails when |
|----------|--------------|-------------|------------|
| Balanced (orthogonal) assignment | every condition appears equally in every batch | always achievable when batches hold >=1 of each condition | not all conditions fit per batch |
| Block randomization across batches | randomize condition within each batch | batch = a block; conditions fit per batch | batch variance is genuinely zero (rare) |
| Incomplete block + batch in model | conditions split across smaller batches, batch term retained | plate/chip smaller than #conditions | unbalanced split inflates artifacts (Nygaard 2016) |
| Reference / bridge sample per batch | shared anchor measured in every batch | cross-batch normalization (TMT proteomics, large cohorts) | anchor not representative |
| Multiplexing + demultiplexing | pool biological units in one lane, split by barcode/genotype | breaking the donor<->lane confound (scRNA-seq) | insufficient SNPs/hashes to assign |
| Run-order randomization | randomize processing/injection order | position/time gradients (LC-MS, plate edge) | order set by convenience |

## Decision Tree by Scenario

| Scenario | Recommended design | Why |
|----------|--------------------|-----|
| 24 samples, 3 batches, 2 conditions | balanced: 4 of each condition per batch | batch orthogonal to condition; estimable |
| Conditions outnumber batch capacity | incomplete block; keep batch in the DE model | preserves estimability; no scrubbing |
| Large cohort across many runs | include a shared reference sample per batch | enables cross-batch normalization |
| scRNA-seq, several donors, few lanes | pool donors per lane, demultiplex (demuxlet / hashing) | removes donor<->lane confound (Kang 2018) |
| Hidden/unknown technical structure suspected | estimate surrogate variables (SVA), include in model | captures unmodeled variation (Leek & Storey 2007) |
| Design already confounds batch with condition | redesign; no correction will rescue it | non-identifiable (Nygaard 2016; Leek 2010) |
| General randomization / blocking / unit choice | -> experimental-design/randomization-blocking | foundational design structure |
| Running ComBat-seq / RUVSeq / SVA on real data | -> differential-expression/batch-correction | execution lives there; this skill decides |

## Confounded vs Balanced -- the Canonical Contrast

**Goal:** Make batch effects correctable by keeping the biological variable orthogonal to batch.

**Approach:** Never place all of one condition in one batch. Distribute conditions (and known covariates such as sex) equally across batches so a linear model can estimate batch and condition separately.

```r
# BAD (confounded): batch is aliased with condition -> non-identifiable
#   batch 1: treat, treat, treat, treat       batch 2: ctrl, ctrl, ctrl, ctrl
# GOOD (balanced): batch is orthogonal to condition -> batch effect estimable, removable
#   batch 1: 2 treat + 2 ctrl                 batch 2: 2 treat + 2 ctrl
```

## Constrained Sample-to-Batch Assignment

**Goal:** Allocate samples to batches/lanes/plates to minimize correlation between batch and the biological variables of interest.

**Approach:** Use a block-randomization-with-optimization tool that scores assignments by how evenly the biological factors spread across batches and returns a near-optimal layout.

```r
library(designit)                              # verify API against installed vignette
samples <- data.frame(id = sprintf('S%02d', 1:24),
                      condition = rep(c('ctrl', 'treat'), each = 12),
                      sex = rep(c('M', 'F'), 12))
bc <- BatchContainer$new(dimensions = list(batch = 3, position = 8))
bc <- assign_in_order(bc, samples = samples)
bc <- optimize_design(
  bc,
  scoring = osat_score_generator(batch_vars = 'batch',
                                 feature_vars = c('condition', 'sex')))   # balance both factors
assignment <- bc$get_samples()                # R6 method on the container (no standalone get_samples())

# OSAT alternative (Bioconductor): build a setup object, then optimal.shuffle() -- NOT a bare osat().
```

## Downstream Correction -- Choose by Design, Execute Elsewhere

Correction method selection is a design decision; the execution lives in differential-expression/batch-correction (and single-cell/batch-integration for scRNA-seq). Prefer keeping batch in the analysis model over producing a "cleaned" matrix for inference.

| Method | When it applies | Assumption / caveat | Owner of execution |
|--------|-----------------|---------------------|--------------------|
| Batch as a model covariate | batch known, balanced | charges df honestly; the default for inference | differential-expression |
| ComBat-seq | known batch, integer counts | batch ~ orthogonal to biology; Nygaard caveat if unbalanced | differential-expression/batch-correction |
| ComBat (parametric eB) | known batch, log-normalized data | Gaussian; not for raw counts | differential-expression/batch-correction |
| RUVSeq (RUVg/RUVs/RUVr) | negative-control genes/samples available | controls must be truly null to the biology | differential-expression/batch-correction |
| SVA | hidden/unknown structure | surrogate variables can absorb biology if confounded | this skill estimates; DE consumes |
| limma removeBatchEffect | visualization/clustering ONLY | not for the hypothesis test | data-visualization |
| Harmony / scVI / Seurat anchors | scRNA-seq integration | integration, not DE inference | single-cell/batch-integration |

## Detecting Hidden Batch Effects (SVA)

**Goal:** Estimate unmodeled technical structure (hidden batches) so it can be included in the downstream model.

**Approach:** Fit a model matrix for the biological variable and a null matrix, estimate the number of surrogate variables, then compute them for inclusion as covariates in the DE analysis.

```r
library(sva)
mod  <- model.matrix(~ condition, data = colData)   # full model
mod0 <- model.matrix(~ 1, data = colData)           # null model
n_sv <- num.sv(expr_normalized, mod)                # estimate number of hidden batches
svobj <- sva(expr_normalized, mod, mod0, n.sv = n_sv)
# Add svobj$sv to the design used by differential-expression/de-results; do NOT subtract them
# from the data for the test (subtracting is for visualization only).
```

## Per-Method Failure Modes

### Batch confounded with condition
- **Trigger:** all of one condition processed in one batch/run.
- **Mechanism:** batch and condition are aliased -> non-identifiable.
- **Symptom:** condition effect vanishes (or an artifact appears) after correction.
- **Fix:** redesign with balanced assignment; no post-hoc method recovers it (Leek 2010).

### ComBat on unbalanced groups
- **Trigger:** ComBat applied when condition is partially confounded with batch.
- **Mechanism:** mean-centering batches injects between-group differences and understates residual variance.
- **Symptom:** inflated DE counts and over-confident downstream inference (Nygaard 2016).
- **Fix:** keep batch in the model (ComBat-seq with the biological covariate, or batch as a DE covariate); best is balanced design.

### Subtracting surrogate variables before testing
- **Trigger:** feeding an SV-"cleaned" matrix into the DE test.
- **Mechanism:** double-counts the adjustment and loses degrees of freedom.
- **Symptom:** anti-conservative p-values.
- **Fix:** include SVs as covariates in the model; reserve cleaned matrices for plots.

### Adjusting for a collider
- **Trigger:** "adjust for everything measured" includes a downstream/common-effect variable.
- **Mechanism:** conditioning on a collider opens a spurious path.
- **Symptom:** associations that appear only after adjustment.
- **Fix:** adjust for confounders (common causes), not mediators or colliders.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Balance every condition equally across batches | Leek 2010 *Nat Rev Genet* 11:733 | makes batch estimable and removable |
| Unbalanced ComBat can inflate DE (>1000 vs 11 in one case) | Nygaard 2016 *Biostatistics* 17:29 | mean-centering injects group differences |
| ~50 SNPs/cell suffice to demultiplex pooled donors | Kang 2018 *Nat Biotechnol* 36:89 | breaks donor<->lane confound |
| Keep batch in the model for inference; clean only for viz | Nygaard 2016; Leek 2010 | honest degrees of freedom |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Condition effect disappears after ComBat | batch confounded with condition | balance at design time |
| Inflated DE list after batch correction | unbalanced ComBat | keep batch in model; ComBat-seq with covariate |
| scRNA-seq donor effect equals lane effect | one donor per lane | pool + demultiplex (demuxlet / hashing) |
| Spurious associations after "adjusting for all metadata" | conditioning on a collider/mediator | adjust only for confounders |
| Cannot reconstruct who/when/which-lot | no metadata captured | record date, lot, operator, lane, position |

## Reproducibility Metadata

Record for every sample, because these become the batch/blocking variables: processing date, reagent and kit lot numbers, operator, instrument/flow-cell/lane and well/plate position, library prep batch, and any protocol deviations. Unrecorded technical variation cannot be modeled or balanced after the fact, and is a leakage source for any downstream machine-learning model (see machine-learning/model-validation).

## References

- Leek JT, Scharpf RB, Bravo HC, Simcha D, Langmead B, Johnson WE, Geman D, Baggerly K, Irizarry RA. 2010. Tackling the widespread and critical impact of batch effects in high-throughput data. *Nat Rev Genet* 11:733-739.
- Nygaard V, Rødland EA, Hovig E. 2016. Methods that remove batch effects while retaining group differences may lead to exaggerated confidence in downstream analyses. *Biostatistics* 17:29-39.
- Johnson WE, Li C, Rabinovic A. 2007. Adjusting batch effects in microarray expression data using empirical Bayes methods. *Biostatistics* 8:118-127.
- Zhang Y, Parmigiani G, Johnson WE. 2020. ComBat-seq: batch effect adjustment for RNA-seq count data. *NAR Genom Bioinform* 2:lqaa078.
- Leek JT, Storey JD. 2007. Capturing heterogeneity in gene expression studies by surrogate variable analysis. *PLoS Genet* 3:e161.
- Gagnon-Bartsch JA, Speed TP. 2012. Using control genes to correct for unwanted variation in microarray data. *Biostatistics* 13:539-552.
- Kang HM, Subramaniam M, Targ S, et al. 2018. Multiplexed droplet single-cell RNA-sequencing using natural genetic variation. *Nat Biotechnol* 36:89-94.
- Yan L, Ma C, Wang D, Hu Q, Qin M, Conroy JM, Sucheston LE, Ambrosone CB, Johnson CS, Wang J, Liu S. 2012. OSAT: a tool for sample-to-batch allocations in genomics experiments. *BMC Genomics* 13:689.

## Related Skills

- randomization-blocking - The experimental unit, randomization, and blocking concepts behind a good batch layout
- power-analysis - Account for blocking/batch factors in the power calculation
- sample-size - Balanced designs assume equal n per group
- multiple-testing - Surrogate variables change the effective number of tests
- differential-expression/batch-correction - Executes ComBat-seq/RUVSeq/SVA on real data
- single-cell/batch-integration - scRNA-seq integration (Harmony, scVI, Seurat anchors)
- machine-learning/model-validation - Batch confounding is a data-leakage source for ML
- clinical-biostatistics/power-and-sample-size - Trial randomization and design (regulated regime)

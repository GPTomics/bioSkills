---
name: bio-differential-expression-batch-correction
description: Handles batch effects in bulk RNA-seq via design-matrix inclusion (the correct path for DE), ComBat/ComBat-seq for visualization, SVA for unknown latent factors, RUVSeq for negative-control-gene-anchored unwanted variation, and limma::removeBatchEffect for plotting only. Encodes the Nygaard 2016 cardinal sin against testing on a batch-corrected matrix, the choice between SVA/RUVg/RUVs/RUVr, the confounding non-identifiability problem, the single-cell boundary (Harmony/MNN are NOT for bulk), and the Goh 2017 harmonization critique. Use when designing a DE analysis with batch structure, troubleshooting batch-dominated PCA, choosing ComBat vs ComBat-seq, handling unknown batch via SVA, integrating across studies, or deciding when (rarely) to subtract batch.
tool_type: r
primary_tool: sva
---

## Version Compatibility

Reference examples tested with: sva 3.50+ (includes ComBat + ComBat_seq), DESeq2 1.42+, edgeR 4.0+, limma 3.58+, RUVSeq 1.36+, ggplot2 3.5+, harmony 1.2+ (single-cell context only)

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Batch Effect Correction

**"Remove the batch effect before DE"** -> Almost always WRONG. Include batch as a covariate in the design formula (`~ batch + condition`) so DESeq2/edgeR/limma model it without subtracting. Subtraction is for visualization only.

## The Single Most Important Modern Insight -- The Nygaard 2016 cardinal sin

Nygaard, Rødland, Hovig 2016 *Biostatistics* 17(1):29-39, "Methods that remove batch effects while retaining group differences may lead to exaggerated confidence in downstream analyses." Translation: **never run ComBat (or ComBat-seq, or `removeBatchEffect`, or SVA-subtract-then-test) and then run DE on the corrected matrix.**

Mechanism: batch-correction methods fit a model `y_ij = alpha + X_ij beta + gamma_i + delta_i epsilon_ij` and subtract the batch terms. The downstream DE test then computes p-values as if those degrees of freedom had never been spent. Residual df is lower than what DESeq2/edgeR/limma assume. Type-I error inflates -- the gene list looks more significant than it should.

The right approach: **include batch in the design**. `~ batch + condition`. DESeq2/edgeR/limma will properly partial out the batch effect from the condition estimate and account for the spent degrees of freedom in the inference. The batch is corrected at the inference stage, not by mutating counts.

`removeBatchEffect` (limma) is for visualization only -- the function's help page says so. ComBat/ComBat-seq output is for visualization, clustering, or downstream tools that cannot take a design matrix (rare; mostly ML).

A second clarification: structural confounding (every "treated" sample in batch 1, every "control" in batch 2) is non-identifiable. No method fixes this -- the "treatment effect" and "batch effect" are mathematically the same vector. The fix is experimental design (randomize batches in advance). ComBat/SVA on a fully-confounded design silently removes the treatment effect along with the batch effect.

## Algorithmic Taxonomy

| Method | Input | Mechanism | Use for |
|--------|-------|-----------|---------|
| Design-matrix inclusion (`~ batch + condition`) | Raw counts; known batch | Partial out batch in the GLM, spent df accounted | DE testing -- the correct path |
| ComBat (Johnson, Li, Rabinovic 2007) | Log-transformed / continuous expression | Empirical-Bayes location and scale shifts per batch | Visualization of microarray / continuous data |
| ComBat-seq (Zhang, Parmigiani, Johnson 2020) | Raw RNA-seq counts | NB-GLM equivalent of ComBat; returns integer counts | Visualization of RNA-seq counts; cross-study harmonization for ML (with caveat) |
| `limma::removeBatchEffect` | Normalized expression | Linear regression subtraction with design protection | Visualization only -- explicit help-page warning |
| SVA (Leek, Storey 2007 + Leek 2012) | Normalized expression | Estimate latent surrogate variables explaining residual variance independent of variable-of-interest | Unknown batch / hidden technical structure -- add SVs to design |
| svaseq (Leek 2014) | Counts | Count-data variant of SVA | Counts version |
| RUVg (Risso, Ngai, Speed, Dudoit 2014) | Counts; negative control genes | Factor analysis of control-gene residual; W as covariate | Strong negative controls (ERCC, housekeeping) -- add W to design |
| RUVs | Counts; replicate samples | Factor analysis within replicate groups; assumes replicates differ only by unwanted variation | Multi-condition with biological replicates as "controls" |
| RUVr | Counts; design only | Factor analysis of residuals from initial fit | Most data-driven; most likely to absorb biology -- caution |
| Harmony, MNN, Scanorama, BBKNN | Single-cell embeddings | Iterative alignment of clusters across samples | Single-cell ONLY; not bulk |

## Decision Tree by Scenario

| Scenario | Recommended approach | Why |
|----------|---------------------|-----|
| Known batch, want DE | `~ batch + condition` in design; do NOT subtract | Cardinal sin avoidance |
| PCA shows batch separation | Include batch in design; for the figure, `removeBatchEffect` is OK (visualization only) | Two purposes, two tools |
| Unknown batch structure | `sva` / `svaseq`; add SVs as covariates in design | Captures latent technical factors |
| Have ERCC spike-ins or trusted housekeeping | `RUVSeq::RUVg` with control gene indices; add W to design | Most principled UV removal |
| Have replicate samples (technical reps within biological) | `RUVSeq::RUVs` | Replicate structure indicates "this differs only by UV" |
| Cross-study integration (TCGA + ICGC + own data) | ComBat-seq for counts, ComBat for log-expression, THEN meta-analysis -- do NOT pool then DE | Goh 2017 warning |
| Fully confounded batch and condition | Re-collect samples; no method fixes this | Non-identifiable |
| Visualizing batch removal for a figure | `removeBatchEffect(expr, batch = batch, design = model.matrix(~condition))` | Visualization is what it's for |
| Single-cell data | Harmony / MNN / Scanorama in the single-cell category, not here | Different problem |

## Standard Workflow -- Design-Matrix Inclusion

**Goal:** Test the condition effect while accounting for known batch structure with proper degrees-of-freedom accounting.

**Approach:** Include batch as a covariate in the design formula; DESeq2/edgeR/limma will partial it out and compute correct p-values.

```r
library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata,
                               design = ~ batch + condition)
dds <- DESeq(dds)
res <- results(dds, name = 'condition_treated_vs_control')
```

```r
library(edgeR)
y <- DGEList(counts = counts, group = coldata$condition)
keep <- filterByExpr(y, design = model.matrix(~ batch + condition, coldata))
y <- y[keep, , keep.lib.sizes = FALSE]
y <- normLibSizes(y)

design <- model.matrix(~ batch + condition, coldata)
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)
qlf <- glmQLFTest(fit, coef = 'conditiontreated')
```

When known confounders are continuous (RIN, library prep date as days), include them as continuous covariates -- no batch correction needed for those:

```r
design = ~ RIN + library_prep_day + condition
```

## ComBat-seq (Visualization or Cross-Study Counts)

**Goal:** Adjust raw counts to remove known batch effects while preserving biology, for VISUALIZATION or for downstream tools that need a single corrected matrix.

**Approach:** `ComBat_seq(counts, batch, group)` returns batch-adjusted integer counts via NB-GLM. Use the output for PCA, clustering, ML -- NOT for DE testing.

```r
library(sva)

corrected_counts <- ComBat_seq(counts = as.matrix(counts),
                                batch = coldata$batch,
                                group = coldata$condition,
                                full_mod = TRUE)

vsd_corrected <- vst(DESeqDataSetFromMatrix(corrected_counts, coldata, ~1))
plotPCA(vsd_corrected, intgroup = 'condition')
```

`full_mod = TRUE` keeps biological covariates protected (the `group` argument). With `full_mod = FALSE`, ComBat-seq removes batch AND group differences -- a serious failure mode.

`ComBat` (the non-seq version) is for log-transformed or microarray data, NOT raw counts. Running `ComBat` on raw counts produces fractional values and assumes Gaussian residuals where counts are NB. Wrong tool, silent failure.

## SVA -- Unknown Batch

**Goal:** Discover latent technical factors when batch is not recorded, add them as covariates.

**Approach:** Estimate the number of surrogate variables; estimate the SVs themselves; add them to the design and re-run DE.

```r
library(sva)
library(DESeq2)

dds <- DESeqDataSetFromMatrix(counts, coldata, design = ~ condition)
dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)

mod  <- model.matrix(~ condition, coldata)
mod0 <- model.matrix(~ 1, coldata)

n_sv <- num.sv(norm_counts, mod, method = 'leek')
svobj <- svaseq(norm_counts, mod, mod0, n.sv = n_sv)

for (i in seq_len(ncol(svobj$sv))) {
    colData(dds)[[paste0('SV', i)]] <- svobj$sv[, i]
}
sv_formula <- as.formula(paste('~', paste(paste0('SV', seq_len(ncol(svobj$sv))),
                                          collapse = ' + '), '+ condition'))
design(dds) <- sv_formula
dds <- DESeq(dds)
```

CRITICAL CHECK: compute the correlation between each SV and the variable of interest. If any SV correlates with treatment at r > 0.3, do NOT include it -- doing so would partial out biology, deflating the effect of interest.

```r
sapply(seq_len(ncol(svobj$sv)),
       function(i) cor(svobj$sv[, i], as.numeric(coldata$condition)))
```

`svaseq` is the count-data version (Leek 2014); `sva` is for log-transformed / microarray data.

## RUVSeq -- Negative-Control-Anchored UV

**Goal:** Estimate unwanted variation using ERCC spike-ins, housekeeping genes, replicate samples, or post-fit residuals; add as covariates in design.

**Approach:** Pick the RUV variant matching the available controls; add the resulting W matrix to the design.

```r
library(RUVSeq)
library(edgeR)

control_idx <- which(rownames(counts) %in% ercc_genes)

set <- newSeqExpressionSet(as.matrix(counts), phenoData = coldata)
set <- betweenLaneNormalization(set, which = 'upper')

ruv <- RUVg(set, control_idx, k = 2)

design <- model.matrix(~ pData(ruv)$W_1 + pData(ruv)$W_2 + condition,
                        data = coldata)
y <- DGEList(counts = counts, group = coldata$condition)
y <- normLibSizes(y)
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)
qlf <- glmQLFTest(fit, coef = 'conditiontreated')
```

| Variant | Control source | Most appropriate when | Failure mode |
|---------|----------------|----------------------|--------------|
| RUVg | Negative control genes (ERCC, housekeeping) | Strong, validated controls exist | Bad controls -> partials out biology |
| RUVs | Replicate samples / technical reps | Multi-condition with reps that should differ only by UV | Wrong replicate structure absorbs condition effect |
| RUVr | Residuals from an initial GLM | No external controls; most data-driven | Most likely to absorb true biology -- use last |

Choice of k (number of unwanted factors): no automatic procedure. Try k = 1, 2, 3 and inspect PCA after RUV correction. Treatment effect should remain visible; if it disappears, k is too high. 5-10 is typical; >10 is suspicious.

## removeBatchEffect -- Visualization Only

```r
library(limma)

design <- model.matrix(~ condition, coldata)
log_expr_corrected <- removeBatchEffect(log_expr,
                                         batch = coldata$batch,
                                         design = design)

library(ggplot2)
pca <- prcomp(t(log_expr_corrected), scale. = TRUE)
ggplot(data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2],
                  condition = coldata$condition, batch = coldata$batch),
       aes(PC1, PC2, color = condition, shape = batch)) +
    geom_point(size = 3) +
    ggtitle('After removeBatchEffect (visualization only)')
```

The `design` argument protects condition while removing batch. Common error: passing batch ALSO in the design ("Coefficients not estimable" warning) -- pass batch via `batch=` only, NOT also in `design=`.

DO NOT feed `log_expr_corrected` back into limma `lmFit` for DE testing. The function's own help page says so. The right approach: `lmFit(log_expr, model.matrix(~ batch + condition, coldata))` -- include batch as a covariate.

## Detecting Confounding Early

```r
ct <- table(coldata$condition, coldata$batch)
ct
```

| Pattern | Status | Action |
|---------|--------|--------|
| All cells > 0 with roughly equal proportions | Balanced | Include batch as covariate |
| Some cells low (e.g., 4/4/2/0) | Partially confounded | Include batch; report reduced power |
| Some cells zero AND one factor entirely in one batch | Perfectly confounded | UNFIXABLE -- re-collect or drop the affected comparison |

`alias()` reveals collinear columns in a design matrix:

```r
alias(model.matrix(~ batch + condition, coldata))$Complete
```

## The Single-Cell Boundary

Harmony (Korsunsky et al. 2019 *Nat Methods* 16:1289), MNN (Haghverdi et al. 2018 *Nat Biotechnol* 36:421), Scanorama, BBKNN are designed for single-cell integration where the goal is aligning cell-type clusters across samples/batches. They modify the embedding (PCA coordinates) for downstream UMAP/clustering -- they assume the structural alignment problem of single-cell (many similar cells; align clusters).

DO NOT apply to bulk RNA-seq. Bulk lacks the cluster structure these methods assume.

Conversely, DO NOT use ComBat / ComBat-seq for single-cell. ComBat assumes much more homogeneity than scRNA exhibits; it does not handle zeros well.

See `single-cell/batch-integration` for the single-cell methods.

## The Goh 2017 Critique on "Harmonization"

Goh, Wang, Wong 2017 *Trends Biotechnol* 35(6):498-507 warn that batch-correction methods can:

- Introduce NEW batch-like artifacts (false positives) when the batch structure is unclear
- Fail when the "most genes not DE" assumption is violated (heat shock, immune activation, viral host shutoff)
- Inflate downstream confidence (echoing Nygaard 2016)

Their recommendation for cross-study work: per-cohort DE first, then meta-analyze the effect sizes (`metafor`, `limma` `treat` across cohorts). DO NOT pool data, harmonize, then re-test.

"Harmonization" in clinical genomics often means ComBat across cohorts then DE on the harmonized matrix -- the exact failure mode Nygaard 2016 describes. Same problem, different vocabulary.

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| Design-inclusion (`~ batch + condition`) and SVA give different gene lists | SVA captured biology (high SV-condition correlation) OR known batch was incomplete | Check SV-vs-condition correlations; if any |r|>0.3, drop that SV. If SVA still differs after pruning, hidden technical factor exists beyond known batch -- combine `~ batch + SV_clean + condition` |
| Design-inclusion and RUVg disagree | RUVg control genes (ERCC, housekeeping) had real biological variation | Validate control set; switch RUVg -> RUVs (replicate-based) or back to design only |
| ComBat-seq corrected gene list larger than design-included gene list | Nygaard 2016 cardinal sin: spent df not accounted for in downstream DE on corrected matrix | DISCARD the ComBat-seq DE result; report design-included only |
| Cross-cohort meta-analysis shows different DE per cohort | Real biological heterogeneity OR cohort-specific technical drift | Per-cohort DE then meta-analyze effect sizes (`metafor`); do NOT pool then test (Goh 2017) |
| All methods (design, SVA, RUVg) agree on top 100 genes | Robust signal | Report the intersection as high-confidence |
| All methods disagree | Confounded design OR weak signal at high noise floor | Suspect non-identifiability; verify `alias()` and balance |

The design-inclusion result (`~ batch + condition` or `~ confounders + condition`) is the reference. SVA/RUV results are sensitivity analyses; when they agree with the design result, confidence rises. When they diverge, investigate before believing either.

## Per-Method Failure Modes

### Tested on a batch-corrected matrix -- inflated significance

**Trigger:** Pipeline: `ComBat_seq(counts, batch)` -> `DESeq2 on corrected_counts`. Many more "significant" genes than expected.

**Mechanism:** ComBat-seq removed batch terms; DESeq2 then computed inference as if those df had not been spent. Type-I error inflates.

**Symptom:** Implausibly long DE gene list; replication in independent cohort recovers <50%; p-value histogram leans anti-conservative.

**Fix:** Re-run with `~ batch + condition` design on RAW counts. Discard the corrected-counts DE result.

### SVA captured the biology

**Trigger:** `sva` returned 5 SVs; SV2 correlates with `condition` at r = 0.7; user added all SVs to design.

**Mechanism:** When biology is correlated with the hidden factor (e.g., disease severity drives both expression AND blood-draw timing), SVA's SVs capture biology along with technical noise. Partialling them out deflates the effect of interest.

**Symptom:** Condition effect that was clear in raw PCA disappears after SV adjustment; few or no DE genes.

**Fix:** Compute SV-vs-condition correlations; exclude any SV with |r| > 0.3 from the design. Re-fit.

### ComBat applied to raw counts

**Trigger:** Pipeline uses `ComBat(counts, batch)` (not `ComBat_seq`); output is fractional.

**Mechanism:** ComBat is for log-transformed / Gaussian data. Counts are NB. The Gaussian assumption is wrong and the output is uninterpretable as counts.

**Symptom:** Corrected matrix has fractional values; DESeq2 errors out ("counts matrix should be integers"); naive `round()` gives garbage.

**Fix:** Use `ComBat_seq()` for counts. Or better, include batch as a design covariate.

### Confounded design corrected with SVA -- biology gone

**Trigger:** All treated samples in batch 1, all control in batch 2; user runs SVA hoping it will rescue the design.

**Mechanism:** Batch and treatment vectors are identical (or nearly so). SVA estimates "the unwanted factor"; "the unwanted factor" is treatment.

**Symptom:** SV1 correlates with treatment at r ~ 1; after adjustment, no DE genes.

**Fix:** Acknowledge the design is non-identifiable. No statistical method fixes structural confounding. Re-collect with randomized batches.

### Used Harmony on bulk RNA-seq

**Trigger:** Bulk RNA-seq with batch effect; user reaches for Harmony because they've used it for single-cell.

**Mechanism:** Harmony aligns cluster centroids in single-cell embeddings. Bulk has no cluster structure -- typically 6-30 samples, not 10000+ cells.

**Symptom:** Harmony "succeeds" but the result is meaningless; downstream DE is nonsense.

**Fix:** Use design-matrix inclusion (`~ batch + condition`) or ComBat-seq for visualization. Harmony belongs to `single-cell/batch-integration`.

## Common errors

| Error / symptom | Cause | Fix |
|-----------------|-------|-----|
| `Coefficients not estimable` from removeBatchEffect | Batch included in both `batch=` and `design=` | Pass batch via `batch=` only |
| ComBat output has fractional values | Wrong tool for counts | Use `ComBat_seq` |
| SV adjustment kills condition effect | SVs captured biology | Check correlation; exclude high-correlation SVs |
| `Inconsistent dimensions` in RUVg | `control_idx` is symbols vs ENSEMBL rownames | Match index type |
| `full_mod` not specified in ComBat-seq | Defaults; biology may not be protected | Set `full_mod = TRUE` and pass `group =` |

## References

- Nygaard V, Rødland EA, Hovig E. 2016. Methods that remove batch effects while retaining group differences may lead to exaggerated confidence in downstream analyses. *Biostatistics* 17(1):29-39. doi:10.1093/biostatistics/kxv027
- Johnson WE, Li C, Rabinovic A. 2007. Adjusting batch effects in microarray expression data using empirical Bayes methods. *Biostatistics* 8(1):118-127. doi:10.1093/biostatistics/kxj037
- Zhang Y, Parmigiani G, Johnson WE. 2020. ComBat-seq: batch effect adjustment for RNA-seq count data. *NAR Genom Bioinform* 2(3):lqaa078. doi:10.1093/nargab/lqaa078
- Leek JT, Storey JD. 2007. Capturing heterogeneity in gene expression studies by surrogate variable analysis. *PLoS Genet* 3(9):e161. doi:10.1371/journal.pgen.0030161
- Leek JT, Johnson WE, Parker HS, Jaffe AE, Storey JD. 2012. The sva package for removing batch effects and other unwanted variation in high-throughput experiments. *Bioinformatics* 28(6):882-883. doi:10.1093/bioinformatics/bts034
- Leek JT. 2014. svaseq: removing batch effects and other unwanted noise from sequencing data. *Nucleic Acids Res* 42(21):e161. doi:10.1093/nar/gku864
- Risso D, Ngai J, Speed TP, Dudoit S. 2014. Normalization of RNA-seq data using factor analysis of control genes or samples. *Nat Biotechnol* 32(9):896-902. doi:10.1038/nbt.2931
- Goh WWB, Wang W, Wong L. 2017. Why batch effects matter in omics data, and how to avoid them. *Trends Biotechnol* 35(6):498-507. doi:10.1016/j.tibtech.2017.02.012
- Korsunsky I et al. 2019. Fast, sensitive and accurate integration of single-cell data with Harmony. *Nat Methods* 16(12):1289-1296. doi:10.1038/s41592-019-0619-0
- Haghverdi L, Lun ATL, Morgan MD, Marioni JC. 2018. Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors. *Nat Biotechnol* 36(5):421-427. doi:10.1038/nbt.4091
- Ritchie ME et al. 2015. limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Res* 43(7):e47. doi:10.1093/nar/gkv007

## Related Skills

- deseq2-basics - DE with `~ batch + condition` design
- edger-basics - edgeR pipeline with batch in design
- de-results - p-value histogram diagnostics for missed batch
- de-visualization - PCA shows batch effects; sample distance heatmap; figure-time use of removeBatchEffect
- expression-matrix/metadata-joins - Confounding detection; sample swap detection
- expression-matrix/normalization - TMM/RLE failure modes overlap with batch failure modes
- single-cell/batch-integration - Harmony, MNN, Scanorama for single-cell (not bulk)

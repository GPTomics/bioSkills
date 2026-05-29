---
name: bio-differential-expression-edger-basics
description: Performs differential expression on bulk RNA-seq count data with edgeR's negative-binomial GLM and quasi-likelihood F-test framework. Covers DGEList construction, filterByExpr, TMM/TMMwsp normalization, robust dispersion estimation, glmQLFit/glmQLFTest, TREAT for magnitude-bounded hypotheses, contrasts via no-intercept designs, voom and voomWithQualityWeights for heterogeneous samples, and the edgeR v4 bias-corrected APL changes. Use when running bulk DE with edgeR, choosing edgeR over DESeq2 (small n, transcript DE via catchSalmon, large samples), needing TREAT for a fold-change-threshold hypothesis, troubleshooting v3-to-v4 reproducibility, building paired or interaction designs, or handling library-quality heterogeneity.
tool_type: r
primary_tool: edgeR
---

## Version Compatibility

Reference examples tested with: edgeR 4.0+, limma 3.58+, statmod 1.5+ (for voom internals), tximport 1.30+ (for catchSalmon path)

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# edgeR Basics

**"Find differentially expressed genes between conditions"** -> Fit a negative-binomial GLM per gene with empirical-Bayes-moderated dispersions, test coefficients with the quasi-likelihood F-test (proper finite-sample FPR control), and report ranked DE lists.

## The Single Most Important Modern Insight -- edgeR v4 changed the QL framework and `legacy=FALSE` is now default

Chen, Chen, Lun, Baldoni, Smyth (2025) *Nucleic Acids Res* 53(2):gkaf018 introduced a bias-corrected adjusted profile likelihood (APL) for dispersion estimation that handles small/zero counts properly, and reworked the QL framework. `glmQLFit()` in v4 takes `legacy=FALSE` by default. Old (v3) pipelines that worked perfectly in 2023 now produce subtly different numbers in 2025 -- not wrong, just different. If a tutorial result doesn't match, check which version was used; set `legacy=TRUE` to reproduce v3 exactly OR (preferred) re-run with v4 defaults and accept the new -- better -- numbers.

Two other v4 changes worth knowing: (1) `calcNormFactors()` is deprecated in favor of `normLibSizes()` (same function, new name); (2) `method='TMM'` remains the documented default; `method='TMMwsp'` (TMM with singleton pairing) is an alternative introduced for samples with many zeros and is the preferred choice for sparse / low-count data. Pass `method=` explicitly for reproducibility. Old names still work, so old scripts run, but new code should use the new names.

The choice between Wald-equivalent (`glmLRT`) and QL-F (`glmQLFTest`) is not a stylistic preference: `glmQLFTest` is the modern default because it accounts for uncertainty in the dispersion estimate via an additional QL dispersion. `glmLRT` is anti-conservative with small n. The two p-values differ -- `glmQLFTest` p-values are typically >= `glmLRT` p-values for the same model (the second-stage QL dispersion correction is usually >= 1).

## Algorithmic Taxonomy

| Test | What it does | When mandatory | Failure mode |
|------|--------------|----------------|--------------|
| QL F-test (`glmQLFit` + `glmQLFTest`) | NB GLM with second-stage QL dispersion to model dispersion uncertainty | DEFAULT for any modern bulk DE | None within design assumptions |
| LRT (`glmFit` + `glmLRT`) | Likelihood-ratio of nested GLMs using point dispersion estimate | Only when n>=10/group AND simple design AND no QL F-test available | Anti-conservative with small n -- inflated false positives |
| Exact test (`exactTest`) | Conditional NB test for two groups, no covariates | Legacy; one-factor two-group ONLY | Cannot adjust for batch or covariates |
| TREAT (`glmTreat`) | Tests H0: |LFC| <= tau vs HA: |LFC| > tau (McCarthy & Smyth 2009) | Want FDR control for "biologically meaningful fold change" claim | Conservative; tau must be pre-specified |
| voom + lmFit + eBayes | Linear model with empirical-Bayes moderation on voom-weighted log-CPM | Heterogeneous library sizes (>3x); samples with widely varying quality | Assumes log-normal-after-weighting; less direct count model |
| voomWithQualityWeights | voom + per-sample quality weights from arrayWeights (Liu 2015) | Some samples markedly worse quality (low RIN, contamination) | With very small n, sample weights are noisy |

## Decision Tree by Scenario

| Scenario | Recommended path | Why |
|----------|------------------|-----|
| Standard bulk RNA-seq, n>=3/group, modern script | `filterByExpr -> normLibSizes -> estimateDisp(robust=TRUE) -> glmQLFit(robust=TRUE) -> glmQLFTest` | Modern default with proper FPR control |
| n=2-3/group | edgeR QL F-test (over DESeq2) | Schurch 2016 *RNA* 22:839 + edgeR v4 changes: QL is the tightest FPR control at small n |
| Library sizes vary >3x or RIN heterogeneous | `voom` or `voomWithQualityWeights` + lmFit + eBayes | Precision weights downweight noisy observations per sample |
| Pre-specified biologically meaningful fold-change threshold | `glmTreat(fit, coef=..., lfc=log2(1.5))` | Post-hoc filtering by |LFC| does NOT control FDR for the magnitude hypothesis |
| Multi-group, "any change" omnibus | `glmQLFTest(fit, coef=2:k)` | Joint F-test on multiple coefficients |
| Multi-group, all pairwise contrasts | `~ 0 + group` parameterization + `makeContrasts` | Clean, readable contrasts |
| Transcript-level DE (DTE) | `catchSalmon()` / `catchKallisto()` -> `DGEList` with overdispersion | Baldoni 2024 *NAR* 52:e13: properly handles inferential variance |
| Cross-tool sanity check | Run DESeq2 in parallel; expect >=70% overlap at top 500 | <60% overlap suggests modeling problem, not tool difference |
| Single-cell pseudobulk | Aggregate counts per donor; standard edgeR QL pipeline | Crowell 2020 *Nat Commun* 11:6077: pseudobulk avoids the FDR inflation of cell-level DE |
| Reproducing pre-2025 result | `glmQLFit(legacy=TRUE)` | edgeR v4 introduced bias-corrected APL with `legacy=FALSE` default; v3 was slightly biased |

## Standard Workflow

**Goal:** Take a raw integer count matrix and group labels to a ranked DE table with proper finite-sample FPR control.

**Approach:** DGEList -> design-aware filtering -> TMM normalization (or TMMwsp for sparse data) -> robust dispersion estimation -> QL fit with robust dispersion shrinkage -> QL F-test on the coefficient of interest -> ranked table.

```r
library(edgeR)
library(limma)

y <- DGEList(counts = counts, group = group)
design <- model.matrix(~ group)

keep <- filterByExpr(y, design)
y <- y[keep, , keep.lib.sizes = FALSE]

y <- normLibSizes(y)
y <- estimateDisp(y, design, robust = TRUE)
plotBCV(y)

fit <- glmQLFit(y, design, robust = TRUE)
qlf <- glmQLFTest(fit, coef = 2)

topTags(qlf, n = 20)
all_de <- topTags(qlf, n = Inf, sort.by = 'none')$table
```

The `robust = TRUE` flags propagate Phipson, Lee, Majewski, Alexander, Smyth (2016) *Ann Appl Stat* 10:946 robust hyperparameter estimation through both the NB dispersion shrinkage (`estimateDisp`) AND the QL dispersion shrinkage (`glmQLFit`). Set BOTH; setting only one is the single most common omission. The default `robust=FALSE` is a relic.

## Filtering -- `filterByExpr` Internals

**Goal:** Remove genes with insufficient expression to reduce noise and multiple-testing burden, in a design-aware way.

**Approach:** `filterByExpr(y, design)` uses the smallest group size from the design matrix to set the minimum number of samples; threshold is CPM >= `min.count / median.lib.size * 1e6` AND total count >= `min.total.count`.

Default parameters: `min.count = 10`, `min.total.count = 15`, `large.n = 10`, `min.prop = 0.7`.

```r
keep <- filterByExpr(y, design)
y <- y[keep, , keep.lib.sizes = FALSE]
```

Filter ONCE, before normalization and dispersion estimation. Filtering after `estimateDisp` invalidates the trended dispersion (the trend was fit on the now-smaller gene set). edgeR has no automatic independent filtering like DESeq2 -- skipping `filterByExpr` is the single biggest reason an edgeR analysis underperforms a comparable DESeq2 analysis.

## Normalization -- TMM/TMMwsp Is an Offset, Not a Division

**Goal:** Correct for library composition bias (the "few genes consume disproportionate reads" problem) via a per-sample scaling factor.

**Approach:** `normLibSizes(y)` defaults to `method='TMM'` in v4 (same as v3); `method='TMMwsp'` is the preferred alternative for sparse / single-cell-pseudobulk data with many zeros. Result is a vector of normalization factors stored in `y$samples$norm.factors`; the effective library size is `lib.size * norm.factors`. Counts are NOT divided.

```r
y <- normLibSizes(y)
y$samples$norm.factors

cpm_vis <- cpm(y, normalized.lib.sizes = TRUE)
log_cpm <- cpm(y, log = TRUE, prior.count = 2)
```

The factor enters the GLM as part of the offset. Running TMM on pre-normalized values is wrong (and silent -- gives wrong numbers without error). `prior.count = 2` is the modern edgeR default for `cpm(log=TRUE)`; smaller priors (0.25) make low-count log values noisy; larger priors (5-10) shrink them toward zero.

When TMM/TMMwsp assumptions fail (see Failure Modes: prokaryotic stress, MYC amplification, viral host shutoff): use `method='upperquartile'` or supply known stable reference genes via offsets.

## QL F-test vs LRT vs Exact Test

```r
fit_ql <- glmQLFit(y, design, robust = TRUE)
qlf <- glmQLFTest(fit_ql, coef = 2)

fit_lrt <- glmFit(y, design)
lrt <- glmLRT(fit_lrt, coef = 2)

y <- estimateDisp(y)
et <- exactTest(y)
```

QL F-test is the modern default. LRT is anti-conservative with small n because it does not account for dispersion uncertainty. The exact test handles only two-group one-factor designs and exists for backward compatibility -- use the GLM pipeline for anything with covariates or blocking.

In v4, the QL F-test is bias-corrected via the new APL when `legacy=FALSE` (the v4 default). To exactly reproduce a v3 result: `glmQLFit(y, design, robust = TRUE, legacy = TRUE)`.

## TREAT -- Testing Against a Fold-Change Threshold

**Goal:** Control FDR for the hypothesis "biologically meaningful fold change", not "fold change non-zero".

**Approach:** `glmTreat(fit, coef=..., lfc=log2(tau))` tests H0: |LFC| <= tau vs HA: |LFC| > tau. The threshold tau MUST be biologically pre-specified -- choosing tau after seeing the data is p-hacking.

```r
tr <- glmTreat(fit, coef = 2, lfc = log2(1.5))
topTags(tr)
```

The cardinal sin TREAT avoids: filtering `padj < 0.05 & abs(logFC) > 1` after a vanilla QL F-test does NOT control FDR for "the LFC exceeds 1". Post-hoc filtering is FDR-controlled for the |LFC|>0 hypothesis only. If reviewers ask "what is the FDR of the '2-fold up' gene list", the only honest answers are TREAT (FDR controlled at the threshold) or "FDR is for non-zero only; the magnitude filter has no FDR guarantee".

Equivalent in DESeq2: `results(dds, lfcThreshold = log2(1.5), altHypothesis = 'greaterAbs')`. Both implement the McCarthy & Smyth 2009 *Bioinformatics* 25:765 idea.

## Contrasts via No-Intercept Designs

**Goal:** Define clean pairwise or multi-group contrasts.

**Approach:** Use `~ 0 + group` parameterization so each coefficient is the group mean; build `makeContrasts(...)` for any pairwise or weighted-average comparison.

```r
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)

con <- makeContrasts(
    TreatedVsControl = treated - control,
    DrugAVsDrugB     = drugA - drugB,
    ATvsBT           = (treated_A - control_A) - (treated_B - control_B),
    levels = design
)

qlf_t_vs_c <- glmQLFTest(fit, contrast = con[, 'TreatedVsControl'])
qlf_interaction <- glmQLFTest(fit, contrast = con[, 'ATvsBT'])
```

`makeContrasts` is more readable than numeric vectors for any non-trivial design.

## voom and voomWithQualityWeights

**Goal:** Use limma's linear-model framework with empirical-Bayes moderation, applying voom precision weights for the mean-variance trend of log-CPM.

**Approach:** voom transforms counts to log-CPM with per-observation precision weights derived from the mean-variance trend; lmFit + eBayes proceeds as for microarrays.

```r
v <- voom(y, design, plot = TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit, robust = TRUE)
tt <- topTable(fit, coef = 2, number = Inf)
```

`voomWithQualityWeights` (Liu, Holik, Su et al. 2015 *NAR* 43:e97) adds per-sample weights to downweight outlier samples. Use when:
- Library sizes vary >5x across samples
- RIN varies >2 units
- PCA shows one or more samples clearly off the main cluster

```r
v <- voomWithQualityWeights(y, design, plot = TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit, robust = TRUE)
```

`eBayes(robust = TRUE)` uses Phipson 2016 robust hyperparameter estimation -- NOT the default but strongly recommended for RNA-seq.

## Transcript-Level DE via catchSalmon

**Goal:** Test differential transcript expression with proper inferential variance from quantification bootstrap replicates.

**Approach:** `catchSalmon()` / `catchKallisto()` import per-transcript counts with overdispersion estimates from the bootstrap replicates; pass to `DGEList`; the rest of the pipeline is standard.

```r
salmon <- catchSalmon(paths = file.path('salmon_out', samples$id))
y_tx <- DGEList(counts = salmon$counts / salmon$annotation$Overdispersion,
                genes  = salmon$annotation)
```

Baldoni, Chen, Hediyeh-zadeh et al. 2024 *NAR* 52:e13 ("Dividing out quantification uncertainty"): the per-transcript Overdispersion column captures the variance contribution from the Salmon EM. Dividing the counts by this scaling provides effective counts that limma/edgeR can model as if quantification was certain.

## edgeR vs DESeq2 vs limma-voom

| Scenario | Recommended | Rationale |
|----------|-------------|-----------|
| Modern default | DESeq2 or edgeR QL (both fine) | ~70-90% overlap on top hits; choose by ecosystem |
| n = 2-3/group | edgeR QL | Tightest FPR at small n |
| Library sizes heterogeneous | limma-voom or voomWithQualityWeights | Precision weights matter most here |
| Salmon/kallisto -> gene-level DGE | DESeq2 (DESeqDataSetFromTximport handles offsets natively) | Cleanest path |
| Salmon -> transcript-level DTE | edgeR `catchSalmon` | Baldoni 2024 framework |
| Need apeglm/ashr LFC shrinkage | DESeq2 | edgeR has no equivalent built-in |
| Need TREAT-style threshold testing | Both (`glmTreat` or `lfcThreshold=`) | Equivalent |
| Python-only environment | PyDESeq2 | No edgeR Python equivalent |

If DESeq2 and edgeR agree on >=70% of the top 500: results are robust. <60%: investigate filtering, normalization, dispersion, or a confounded covariate -- usually a modeling issue.

## Per-Method Failure Modes

### Forgot `filterByExpr` -- inflated multiple-testing burden

**Trigger:** edgeR pipeline run on a count matrix with all genes (~60k for human) including many with zero or near-zero counts; significant gene list smaller than expected.

**Mechanism:** edgeR's QL pipeline does NOT have DESeq2's automatic independent filtering. Every gene tested contributes to the BH denominator.

**Symptom:** Low-count genes dominate the tested set; tested-gene count is 60k instead of ~15-20k; padj distribution skewed toward 1.

**Fix:** `filterByExpr(y, design)` BEFORE normalization and dispersion estimation. Filtering after `estimateDisp` invalidates the trended dispersion.

### Used `glmLRT` with n=3 -- inflated false positives

**Trigger:** Tutorial copy-paste of `glmFit` + `glmLRT` on small-n data; many "significant" genes that don't replicate.

**Mechanism:** `glmLRT` uses a point dispersion estimate; with small n, dispersion is uncertain and the test is anti-conservative.

**Symptom:** Many more DE genes than DESeq2 or limma-voom on the same data; p-value histogram has anti-conservative left tail.

**Fix:** Use `glmQLFit` + `glmQLFTest` with `robust = TRUE` on both.

### TMM/RLE assumption broken -- size factors absorb biology

**Trigger:** Bacterial stress response, viral infection with host shutoff, MYC-amplified tumor vs normal -- biological systems where >50% of genes truly change.

**Mechanism:** TMM/TMMwsp assumes most genes are unchanged. When violated, the "trimmed reference" is dominated by DE genes; the size factor compensates by absorbing the biological shift.

**Symptom:** MA plot shows the bulk cloud shifted off zero; reported fold changes don't match orthogonal validation (qPCR, Western); known DE genes show muted LFC.

**Fix:** `normLibSizes(y, method='upperquartile')` is a partial fix; supply ERCC spike-in offsets or curated stable housekeeping genes via `y$offset` for the principled solution.

### Reproducibility break between edgeR v3 and v4

**Trigger:** Pre-2025 script run on edgeR 4.0+ produces different DE genes than the published result.

**Mechanism:** v4 changed the QL framework (bias-corrected APL) and the default behavior of `glmQLFit` (`legacy = FALSE`). TMM remains the documented default for `normLibSizes`; TMMwsp is an added alternative for sparse data.

**Symptom:** Different significant gene set; absolute LFC differences of 5-15% on low-count genes; reviewer confusion.

**Fix:** For exact reproduction: `glmQLFit(y, design, robust=TRUE, legacy=TRUE)`. The normalization default has not actually changed -- TMM remains the v4 default. Better: rerun and accept the v4 -- improved -- numbers, noting the version change in methods.

## Common errors

| Error / symptom | Cause | Fix |
|-----------------|-------|-----|
| `decidetestsDGE` not found | Removed in edgeR v4 | Use `decideTests(qlf)` |
| `design matrix not full rank` | Confounded covariates | Inspect with `alias(design)$Complete` |
| `No residual df` | Too few replicates for the model | Reduce model complexity or get more samples |
| Script expects `$adj.P.Val` but `topTags` returns `$FDR` | Tool-name column mix-up | edgeR uses `$FDR`; limma uses `$adj.P.Val`; DESeq2 uses `$padj` |
| QL F p-values numerically different from a 2023 paper | edgeR v4 default `legacy=FALSE` | Set `legacy=TRUE` to reproduce v3, or note the version change |

## References

- Robinson MD, McCarthy DJ, Smyth GK. 2010. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. *Bioinformatics* 26(1):139-140. doi:10.1093/bioinformatics/btp616
- McCarthy DJ, Chen Y, Smyth GK. 2012. Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation. *Nucleic Acids Res* 40(10):4288-4297. doi:10.1093/nar/gks042
- Chen Y, Lun ATL, Smyth GK. 2016. From reads to genes to pathways: differential expression analysis of RNA-Seq experiments using Rsubread and the edgeR quasi-likelihood pipeline. *F1000Res* 5:1438. doi:10.12688/f1000research.8987.2
- Chen Y, Chen L, Lun ATL, Baldoni PL, Smyth GK. 2025. edgeR v4: powerful differential analysis of sequencing data with expanded functionality and improved support for small counts and larger datasets. *Nucleic Acids Res* 53(2):gkaf018. doi:10.1093/nar/gkaf018
- Lund SP, Nettleton D, McCarthy DJ, Smyth GK. 2012. Detecting differential expression in RNA-sequence data using quasi-likelihood with shrunken dispersion estimates. *Stat Appl Genet Mol Biol* 11(5):Article 8. doi:10.1515/1544-6115.1826
- Phipson B, Lee S, Majewski IJ, Alexander WS, Smyth GK. 2016. Robust hyperparameter estimation protects against hypervariable genes and improves power to detect differential expression. *Ann Appl Stat* 10(2):946-963. doi:10.1214/16-AOAS920
- Law CW, Chen Y, Shi W, Smyth GK. 2014. voom: precision weights unlock linear model analysis tools for RNA-seq read counts. *Genome Biol* 15(2):R29. doi:10.1186/gb-2014-15-2-r29
- Liu R, Holik AZ, Su S, Jansz N, Chen K, Leong HS, Blewitt ME, Asselin-Labat M-L, Smyth GK, Ritchie ME. 2015. Why weight? Modelling sample and observational level variability improves power in RNA-seq analyses. *Nucleic Acids Res* 43(15):e97. doi:10.1093/nar/gkv412
- McCarthy DJ, Smyth GK. 2009. Testing significance relative to a fold-change threshold is a TREAT. *Bioinformatics* 25(6):765-771. doi:10.1093/bioinformatics/btp053
- Baldoni PL, Chen Y, Hediyeh-zadeh S, Liao Y, Dong X, Ritchie ME, Shi W, Smyth GK. 2024. Dividing out quantification uncertainty allows efficient assessment of differential transcript expression with edgeR. *Nucleic Acids Res* 52(3):e13. doi:10.1093/nar/gkad1167
- Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK. 2015. limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Res* 43(7):e47. doi:10.1093/nar/gkv007
- Schurch NJ et al. 2016. How many biological replicates are needed in an RNA-seq experiment and which differential expression tool should you use? *RNA* 22(6):839-851. doi:10.1261/rna.053959.115
- Robinson MD, Oshlack A. 2010. A scaling normalization method for differential expression analysis of RNA-seq data. *Genome Biol* 11(3):R25. doi:10.1186/gb-2010-11-3-r25

## Related Skills

- deseq2-basics - Cross-check, or use when LFC shrinkage (apeglm) needed
- de-results - FDR, IHW, GSEA preparation, padj column name variation
- de-visualization - BCV plot, MD plot, PCA via plotMDS
- batch-correction - Include batch in design vs correct-then-test cardinal sin
- timeseries-de - voom + splines for time-course
- expression-matrix/counts-ingest - Salmon/kallisto input via tximport or catchSalmon
- expression-matrix/normalization - TMM/TMMwsp/RLE mechanics
- expression-matrix/metadata-joins - Reference level, paired design, interaction parameterization
- expression-matrix/gene-id-mapping - Annotating DE results with symbols
- rna-quantification/tximport-workflow - Detailed tximport mechanics
- pathway-analysis/gsea - Ranked-list input from DE

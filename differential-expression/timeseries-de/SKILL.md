---
name: bio-differential-expression-timeseries-de
description: Analyzes time-series and longitudinal RNA-seq for differential expression and trajectory structure. Covers DESeq2 LRT with reduced models, time as factor vs continuous vs natural splines, maSigPro (Nueda 2014 for RNA-seq), ImpulseDE2 with explicit impulse-model failure modes, DREAM for repeated measures via linear mixed models, pseudoreplication avoidance, conditional vs marginal modeling, and trajectory clustering with DPGP, Mfuzz (with Schwämmle 2010 fuzzifier estimation), and splines+k-means. Use when modeling time-course or longitudinal expression, choosing factor vs spline, handling repeated measures from the same subject, avoiding pseudoreplication, clustering temporal trajectories, or selecting between dedicated time-course tools and pairwise+LRT.
tool_type: r
primary_tool: DESeq2
---

## Version Compatibility

Reference examples tested with: DESeq2 1.42+, edgeR 4.0+, limma 3.58+, splines (base R), maSigPro 1.74+, ImpulseDE2 1.10+ (Bioconductor archive; verify availability), variancePartition / dream 1.32+, Mfuzz 2.62+, TCseq 1.26+, ggplot2 3.5+, pheatmap 1.0+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Time-Series Differential Expression

**"Find genes that change over time"** -> Define what "change" means -- any non-zero time effect (LRT), a smooth nonlinear trend (splines), a transient impulse (ImpulseDE2), or differing trajectories between groups (interaction LRT) -- and choose the model that asks that specific question while handling repeated-measures correctly.

## The Single Most Important Modern Insight -- Most dedicated time-course tools UNDERPERFORM pairwise + LRT on short series

Spies, Renz, Beyer, Ciaudo 2019 *Brief Bioinform* 20:288 benchmarked dedicated time-course tools (ImpulseDE2, splineTC, maSigPro, EBSeqHMM, TimeReg) against naive DESeq2/edgeR pairwise comparisons + LRT for omnibus, on simulated and real time-courses. Finding: **on short series (<8 time points), naive pairwise pattern-of-significance OUTPERFORMS dedicated TC tools** because of high false-positive rates in the latter. The exception is ImpulseDE2, which holds up better than the others -- IF its impulse assumption (rise-then-plateau or fall-then-plateau) actually fits the biology.

For most experimental time courses (3-6 time points, common in pharmacology and developmental biology), the right tool is DESeq2 with `test='LRT'` and a sensible reduced model. Reserve splines for >5 evenly-spaced time points; reserve ImpulseDE2 for monotonic-then-saturating dynamics; reserve DREAM for repeated measures.

A second insight that is constantly violated: pseudoreplication. If 3 subjects each contribute 4 time points (12 samples), the effective sample size for testing TIME effects is closer to 3, not 12 -- the within-subject observations are not independent. Treating them as independent inflates type-I error dramatically. Either include subject as a fixed effect, use DREAM (mixed model), or collapse to per-subject means (loses time info).

## Algorithmic Taxonomy

| Method | What it tests | Best for | Failure mode |
|--------|---------------|----------|--------------|
| DESeq2 LRT (`test='LRT'`, `reduced=`) | Joint hypothesis: dropped terms are jointly zero | Default for "any time effect"; multi-group time interaction | Reports LFC of last coefficient, not omnibus -- use padj only |
| DESeq2 + splines (`ns(time, df=3)`) | Smooth nonlinear time effect | 5+ time points, smooth dynamics, multi-group interaction | Spline df > unique time points / 2 overfits |
| maSigPro (Nueda 2014 RNA-seq update) | Polynomial regression on time per group | Multi-group time-course, regression-style hypotheses | Polynomial assumption can be wrong; less popular than DESeq2 LRT |
| ImpulseDE2 (Fischer, Theis, Yosef 2018) | Constant vs monotonic vs impulse trajectories | Monotonic-then-saturating or impulse-like responses | Fails on oscillatory, multi-phase, monotonic-non-asymptotic; sensitive to noise on short series |
| DREAM (Hoffman, Roussos 2021) | Per-gene linear mixed model with random subject | Repeated measures with >2 time points per subject | Slower; requires variancePartition stack; `ddf='adaptive'` default (Kenward-Roger for n<=20, Satterthwaite otherwise) |
| voom + duplicateCorrelation | Single average within-subject correlation | Technical reps within bio reps, paired pre/post | Single correlation across all genes is approximation |
| edgeR LRT with subject as factor | Subject as fixed effect | Small subject count, simple design | Wastes df; can't handle continuous time well |
| TCseq | Spline-based DE + fuzzy clustering | Combined pipeline for DE + cluster discovery | Single-tool dependency; verify maintenance |

## Decision Tree by Scenario

| Scenario | Recommended approach | Why |
|----------|---------------------|-----|
| 2-3 discrete time points, independent samples per time | DESeq2 LRT, time as factor, `reduced = ~1` | Too few points for splines |
| 5+ time points, single group, smooth dynamics | DESeq2 LRT with `ns(time, df=3)`, `reduced = ~1` | Splines capture nonlinearity efficiently |
| 5+ time points, two groups, "do trajectories differ?" | LRT with `~ treatment * ns(time, df=3)` vs `reduced = ~ treatment + ns(time, df=3)` | Tests the interaction (treatment-specific time response) |
| Same subject sampled over time (longitudinal) | DREAM (random subject) OR DESeq2 with subject fixed effect | Mandatory: pseudoreplication otherwise |
| Monotonic-then-saturating biology expected | ImpulseDE2 | Built for the assumption |
| Circadian or cell-cycle (cyclical) | Fourier basis (`fda::create.fourier.basis`) or dedicated tools (JTK_CYCLE, MetaCycle) | Splines can't represent periodicity |
| Short series (<8 time points) | DESeq2 pairwise + LRT (Spies 2019 finding) | Dedicated TC tools have high FPR on short series |
| Want to cluster trajectories after DE | DPGP (nonparametric), Mfuzz (with Schwämmle 2010 m estimation), splines + k-means | Standardize per gene first |
| Multi-batch time course | Add batch to design; test the interaction term | Standard DESeq2 / edgeR pattern |

## DESeq2 LRT for Time -- The Canonical Pattern

**Goal:** Test whether time has ANY effect (omnibus), or whether time trajectories differ between groups (interaction).

**Approach:** Specify a full design including the time terms; specify a reduced design dropping those terms; LRT compares.

```r
library(DESeq2)

dds <- DESeqDataSetFromMatrix(counts, colData, design = ~ time)
dds <- DESeq(dds, test = 'LRT', reduced = ~ 1)
res <- results(dds)
```

The LRT p-value tests "any difference among time points". The `log2FoldChange` column reports the LAST coefficient in `resultsNames(dds)` -- NOT the omnibus effect. Use padj only from LRT results; extract individual Wald per time point for effect sizes.

Interaction (treatment-specific time response):

```r
dds <- DESeqDataSetFromMatrix(counts, colData,
    design = ~ treatment + time + treatment:time)
dds <- DESeq(dds, test = 'LRT', reduced = ~ treatment + time)
res_interaction <- results(dds)
```

This tests "does the time trajectory differ between treatment groups?" -- a different question than "is there a time effect" or "is there a treatment effect".

## Time as Factor vs Continuous vs Spline

| Encoding | Assumption | df spent | When |
|----------|------------|----------|------|
| Factor (`time` as factor) | No structure; each level independent | (n_levels - 1) | Few discrete time points, irregular spacing, interest in specific pairwise comparisons |
| Continuous (`as.numeric(time)`) | Linear effect on log expression | 1 | Linear biology, well-spaced time points -- often wrong |
| Natural spline (`ns(time, df=k)`) | Smooth nonlinear; k basis functions | k | Dense time courses (5+ points), smooth biology |

Rule of thumb: `df <= unique_time_points / 2`.

```r
library(splines)

dds <- DESeqDataSetFromMatrix(counts, colData,
    design = ~ treatment * ns(time, df = 3))
dds <- DESeq(dds, test = 'LRT', reduced = ~ treatment + ns(time, df = 3))
```

## limma-voom + Splines (Modern Alternative)

```r
library(limma)
library(edgeR)
library(splines)

y <- DGEList(counts = counts)
y <- normLibSizes(y)
keep <- filterByExpr(y, group = metadata$treatment)
y <- y[keep, , keep.lib.sizes = FALSE]

design <- model.matrix(~ treatment * ns(time, df = 3), data = metadata)
v <- voom(y, design, plot = TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit, robust = TRUE)

interaction_cols <- grep(':ns\\(time', colnames(design))
tt <- topTable(fit, coef = interaction_cols, number = Inf)
```

Tests the joint significance of all interaction spline coefficients.

## DREAM for Repeated Measures

**Goal:** Properly model longitudinal data where the same subject is sampled multiple times.

**Approach:** Linear mixed model per gene with subject as a random intercept (or slope), via `variancePartition::dream`. Uses voom-weighted linear regression internally. Default `ddf='adaptive'` uses Kenward-Roger for n <= 20 samples (the typical longitudinal-DE regime) and Satterthwaite otherwise; force either explicitly with `ddf='Kenward-Roger'` or `ddf='Satterthwaite'`.

```r
library(variancePartition)
library(edgeR)
library(BiocParallel)

y <- DGEList(counts = counts)
y <- normLibSizes(y)
keep <- filterByExpr(y, group = metadata$treatment)
y <- y[keep, , keep.lib.sizes = FALSE]

formula <- ~ treatment + time + treatment:time + (1 | subject)

vobj <- voomWithDreamWeights(y, formula, metadata)
fitmm <- dream(vobj, formula, metadata)
fitmm <- eBayes(fitmm)
tt <- topTable(fitmm, coef = grep(':time', colnames(coefficients(fitmm))),
               number = Inf)
```

When to use DREAM over `~ subject + treatment + time` (subject as fixed effect):
- More than 2 time points per subject (more random-effect signal to estimate)
- Many subjects (random effects more parsimonious than fixed)
- Random slopes needed (e.g., `(1 + time | subject)`)
- Multiple random effects (e.g., `(1 | donor) + (1 | batch)`)

`duplicateCorrelation` (`limma`) is the older approximate alternative -- assumes a SINGLE within-subject correlation across all genes. Adequate for technical replicates within biological replicates; less so for proper longitudinal data with multiple time points.

## maSigPro (Nueda 2014 for RNA-seq)

```r
library(maSigPro)

edesign <- data.frame(
    Time      = metadata$time,
    Replicate = metadata$replicate,
    Control   = as.numeric(metadata$treatment == 'Control'),
    Treatment = as.numeric(metadata$treatment == 'Treatment')
)
rownames(edesign) <- metadata$sample

design <- make.design.matrix(edesign, degree = 3)

fit <- p.vector(norm_counts, design, Q = 0.05, MT.adjust = 'BH')
tstep <- T.fit(fit, step.method = 'backward', alfa = 0.05)
sigs <- get.siggenes(tstep, rsq = 0.6, vars = 'groups')
see.genes(sigs$sig.genes, show.fit = TRUE, dis = design$dis,
          cluster.method = 'hclust', k = 9)
```

CITATION: Nueda MJ, Tarazona S, Conesa A (2014) "Next maSigPro: updating maSigPro bioconductor package for RNA-seq time series." *Bioinformatics* 30(18):2598-2602. The earlier Conesa et al. 2006 *Bioinformatics* 22:1096 is the microarray-era original. Cite 2014 for RNA-seq use; many secondary refs miscite 2006.

## ImpulseDE2

**Goal:** Detect transient impulse-like expression patterns (rise then decay, or constant-then-monotonic).

**Approach:** Fits constant, monotonic, and impulse (6-parameter sigmoid: baseline, peak, post-peak baseline, rise rate, decay rate) models per gene; selects the best-fitting and tests for differential dynamics.

```r
library(ImpulseDE2)

dfAnnotation <- data.frame(
    Sample    = colnames(counts),
    Time      = metadata$time,
    Condition = metadata$condition,
    Batch     = metadata$batch
)

imp <- runImpulseDE2(
    matCountData    = as.matrix(counts),
    dfAnnotation    = dfAnnotation,
    boolCaseCtrl    = TRUE,
    vecConfounders  = c('Batch'),
    scaNProc        = 4
)

sig <- imp$dfImpulseDE2Results[imp$dfImpulseDE2Results$padj < 0.05, ]
```

The impulse model FAILS on:
- Oscillatory expression (circadian, cell cycle) -- model can't represent periodicity
- Monotonic-but-non-asymptotic responses (linear increases that don't plateau)
- Multi-phase responses (rise-fall-rise)

For oscillatory biology, use Fourier basis or dedicated tools (JTK_CYCLE, MetaCycle). For non-asymptotic monotonic, splines fit better.

NOTE: ImpulseDE2 was removed from Bioconductor at the 3.13 release (May 2021); last hosted version was 3.10 -- install from the BiocArchive (pin to Bioc 3.10) or the YosefLab GitHub mirror. The Spies 2019 benchmark showed ImpulseDE2 holds up on longer time courses (~8+ time points) but is sensitive to noise on short series; below ~8 time points, DESeq2 pairwise + LRT typically outperforms.

## Trajectory Clustering of DE Genes

**Goal:** After identifying DE-over-time genes, group them by trajectory shape.

**Approach:** Standardize per gene (subtract mean, divide by SD per gene) -- otherwise clusters reflect mean level, not shape. Then apply DPGP (nonparametric), Mfuzz (fuzzy c-means), or splines + k-means.

```r
library(Mfuzz)

eset <- ExpressionSet(assayData = as.matrix(norm_counts[sig_genes, ]))
eset_std <- standardise(eset)

m <- mestimate(eset_std)
cl <- mfuzz(eset_std, c = 9, m = m)

mfuzz.plot(eset_std, cl, mfrow = c(3, 3))
```

The Mfuzz fuzzifier `m` is critical -- too low gives crisp clusters (loses fuzzy advantage); too high collapses everything. `mestimate()` implements Schwämmle & Jensen 2010 *Bioinformatics* 26:2841 to estimate `m` from the data.

CITATION CARE: the Mfuzz PACKAGE paper is Kumar L, Futschik ME (2007) *Bioinformation* 2(1):5-7 (the journal is *Bioinformation*, NOT *Bioinformatics*). The Schwämmle 2010 paper is the fuzzifier-estimation methodology paper, *Bioinformatics*. Many references confuse the two.

For DPGP (Dirichlet Process Gaussian Process; nonparametric in cluster number AND trajectory shape):

CITATION: McDowell IC, Manandhar D, Vockley CM, Schmid AK, Reddy TE, Engelhardt BE (2018) "Clustering gene expression time series data using an infinite Gaussian process mixture model." *PLoS Comput Biol* 14(1):e1005896. CITATION CARE: this is *PLoS Comp Biol*, NOT *Genome Research* (a common miscitation).

Splines + k-means (fast alternative):

```r
library(splines)
spline_coefs <- t(apply(norm_counts[sig_genes, ], 1, function(x) {
    fit <- lm(x ~ ns(metadata$time, df = 4))
    coef(fit)
}))
km <- kmeans(scale(spline_coefs), centers = 6)
```

## Time as the Only Variable (Developmental Series)

When time is the sole variable (e.g., embryonic development series), the DE question becomes "which genes change across the trajectory":

```r
dds <- DESeqDataSetFromMatrix(counts, colData, design = ~ time)
dds <- DESeq(dds, test = 'LRT', reduced = ~ 1)
res <- results(dds)
```

Caveats:
- Reference time point biases the reported LFC (it's the LFC at the LAST time point vs reference)
- For cyclical biology (circadian, cell cycle), use periodic basis functions

## Per-Method Failure Modes

### Pseudoreplication -- treated 12 samples from 3 subjects as 12 independent

**Trigger:** 3 subjects x 4 time points = 12 samples; vanilla DESeq2 with `~ time`; many DE genes.

**Mechanism:** Within-subject observations are correlated; treating them as independent inflates effective sample size from 3 to 12 in the inference, deflating standard errors.

**Symptom:** p-value histogram anti-conservative; many false-positive DE genes; replication fails.

**Fix:** Include subject in design (`~ subject + time`) OR use DREAM with random subject. Collapsing to per-subject means is OK but loses time info.

### LRT reports the wrong LFC

**Trigger:** `DESeq(dds, test='LRT', reduced=~1)` on a 5-time-point study; user reports the `log2FoldChange` column as "the time effect".

**Mechanism:** LRT padj is the omnibus joint test. The LFC reported is for the LAST coefficient in `resultsNames(dds)`, one specific level-vs-reference comparison.

**Symptom:** A 5-time-point factor produces one LFC per gene; reviewer asks "the effect of which time?"

**Fix:** Treat LRT padj as a screen for "any change". For effect sizes, extract Wald coefficients per time point via `results(dds, name='time_T2_vs_T0')` etc.

### ImpulseDE2 reports noise as "impulse"

**Trigger:** Short series (4-5 time points), oscillatory or non-monotonic biology; ImpulseDE2 flags many "impulse" genes that look like noise on inspection.

**Mechanism:** Impulse model has 6 parameters; on short series, easy to fit by chance. For oscillatory data, model is wrong.

**Symptom:** Validation orthogonal data shows the "impulse" genes are not actually transient; replication low.

**Fix:** Use DESeq2 LRT + spline interaction for short or non-impulse data. Reserve ImpulseDE2 for cases where biology is known to be monotonic-then-asymptotic (immune response, cytokine release).

### Wrong spline df

**Trigger:** 4 time points; user sets `ns(time, df = 5)`; bizarre fits per gene.

**Mechanism:** df > number of unique time points produces overfitting; spline basis is rank-deficient.

**Symptom:** Errors about singular fits; or apparently-clean fits that don't generalize.

**Fix:** `df <= unique_time_points / 2`. For 4 time points, df = 2 (or use factor encoding instead).

### Trajectory clusters dominated by expression level

**Trigger:** Mfuzz / k-means clustering of trajectories; clusters separate high- vs low-expressed genes rather than shape patterns.

**Mechanism:** Forgot to standardize per gene; absolute levels dominate distance computations.

**Symptom:** Clusters labeled by mean expression, not trajectory shape.

**Fix:** `standardise()` in Mfuzz, or `scale()` per gene before k-means. The point of trajectory clustering is shape, not magnitude.

### maSigPro miscited as 2006

**Trigger:** Methods section cites "Conesa 2006 maSigPro" for an RNA-seq analysis.

**Mechanism:** The 2006 Conesa paper is the original microarray maSigPro. The 2014 Nueda paper is the RNA-seq update with NB GLM.

**Symptom:** Reviewer asks for the RNA-seq citation specifically.

**Fix:** Cite Nueda MJ, Tarazona S, Conesa A (2014) *Bioinformatics* 30(18):2598-2602 for RNA-seq use.

## Common errors

| Error / symptom | Cause | Fix |
|-----------------|-------|-----|
| `singular fit` from spline model | df > unique time points | Reduce df or use factor encoding |
| LRT p-values numerically identical across genes | Reduced model matches full model (no effect tested) | Verify `reduced` actually drops the term of interest |
| ImpulseDE2 not installable from Bioconductor | Removed at Bioconductor 3.13 (May 2021); last hosted version 3.10 | Install from BiocArchive (pin Bioc 3.10) or YosefLab GitHub mirror |
| Mfuzz clusters dominated by mean level | Forgot standardisation | Use `standardise()` before `mfuzz()` |
| Trajectory plot shows flat lines | Counts not log-transformed before clustering | Use `cpm(y, log=TRUE)` or `vst()` then standardize |
| DREAM very slow | Large gene set with mixed model per gene | Filter to DE-over-time first (LRT screen), then DREAM on the subset |

## References

- Love MI, Huber W, Anders S. 2014. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol* 15(12):550. doi:10.1186/s13059-014-0550-8
- Nueda MJ, Tarazona S, Conesa A. 2014. Next maSigPro: updating maSigPro bioconductor package for RNA-seq time series. *Bioinformatics* 30(18):2598-2602. doi:10.1093/bioinformatics/btu333
- Conesa A, Nueda MJ, Ferrer A, Talón M. 2006. maSigPro: a method to identify significantly differential expression profiles in time-course microarray experiments. *Bioinformatics* 22(9):1096-1102. doi:10.1093/bioinformatics/btl056
- Fischer DS, Theis FJ, Yosef N. 2018. Impulse model-based differential expression analysis of time course sequencing data. *Nucleic Acids Res* 46(20):e119. doi:10.1093/nar/gky675
- Hoffman GE, Roussos P. 2021. dream: powerful differential expression analysis for repeated measures designs. *Bioinformatics* 37(2):192-201. doi:10.1093/bioinformatics/btaa687
- Hoffman GE, Schadt EE. 2016. variancePartition: interpreting drivers of variation in complex gene expression studies. *BMC Bioinformatics* 17:483. doi:10.1186/s12859-016-1323-z
- Spies D, Renz PF, Beyer TA, Ciaudo C. 2019. Comparative analysis of differential gene expression tools for RNA sequencing time course data. *Brief Bioinform* 20(1):288-298. doi:10.1093/bib/bbx115
- Kumar L, Futschik ME. 2007. Mfuzz: a software package for soft clustering of microarray data. *Bioinformation* 2(1):5-7.
- Schwämmle V, Jensen ON. 2010. A simple and fast method to determine the parameters for fuzzy c-means cluster analysis. *Bioinformatics* 26(22):2841-2848. doi:10.1093/bioinformatics/btq534
- McDowell IC, Manandhar D, Vockley CM, Schmid AK, Reddy TE, Engelhardt BE. 2018. Clustering gene expression time series data using an infinite Gaussian process mixture model. *PLoS Comput Biol* 14(1):e1005896. doi:10.1371/journal.pcbi.1005896
- Law CW, Chen Y, Shi W, Smyth GK. 2014. voom: precision weights unlock linear model analysis tools for RNA-seq read counts. *Genome Biol* 15(2):R29. doi:10.1186/gb-2014-15-2-r29

## Related Skills

- deseq2-basics - LRT mechanics; design formulas
- edger-basics - voom + splines pattern
- de-results - LRT padj interpretation; pseudoreplication detection
- de-visualization - Per-gene trajectories; heatmap with time order
- batch-correction - Time-batch confounding; multi-batch time courses
- expression-matrix/metadata-joins - Subject as covariate; repeated measures designs
- pathway-analysis/go-enrichment - Functional analysis of trajectory clusters
- temporal-genomics/circadian-rhythms - Circadian-specific detection (JTK_CYCLE, MetaCycle)
- temporal-genomics/temporal-clustering - Standalone trajectory clustering methods
- temporal-genomics/trajectory-modeling - GAM trajectory fitting
- temporal-genomics/temporal-grn - Dynamic gene regulatory network inference

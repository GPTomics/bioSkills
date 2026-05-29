---
name: bio-differential-expression-deseq2-basics
description: Performs differential expression on bulk RNA-seq count data with DESeq2's negative-binomial GLM, Wald and LRT testing, apeglm/ashr/normal LFC shrinkage, independent filtering, Cook's outlier handling, VST/rlog transforms, and design formulas including paired, batch, and interaction terms. Use when running bulk DE, choosing DESeq2 over edgeR or limma-voom, building a paired or interaction design, applying LFC shrinkage for ranking or GSEA, choosing Wald vs LRT, troubleshooting padj=NA, picking VST vs rlog, importing salmon/kallisto via tximport, or analyzing prokaryotic RNA-seq.
tool_type: r
primary_tool: DESeq2
---

## Version Compatibility

Reference examples tested with: DESeq2 1.42+, apeglm 1.28+, ashr 2.2+, IHW 1.34+, tximport 1.30+, edgeR 4.0+ (for cross-comparison), PyDESeq2 0.5+

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# DESeq2 Basics

**"Find differentially expressed genes between conditions"** -> Fit a negative-binomial GLM per gene with shared dispersion shrinkage, test the coefficient of interest (Wald) or the joint effect of a factor (LRT), and report a shrunken effect-size estimate for ranking.

## The Single Most Important Modern Insight -- Shrunken LFC and the Wald p-value come from different models

`lfcShrink()` returns LFCs from a Bayesian posterior with apeglm/ashr/normal priors, BUT the p-value column it carries forward is still the **unshrunken Wald p-value** from `results()`. This is a deliberate design choice (Zhu, Ibrahim, Love 2019 *Bioinformatics* 35:2084) -- the shrunken estimate is for ranking and visualization; the p-value is for inference. Reporting "shrunken LFC = 0.4, padj = 1e-8" mixes two models, which is fine because both are correct for their stated purpose. What is NOT fine: using the shrunken LFC in a downstream filter and then claiming FDR control on that filter (it has none). For threshold-based FDR claims, use `lfcThreshold=` or TREAT (`glmTreat` in edgeR).

A second consequence: `results(dds)` with no `name=` or `contrast=` argument silently returns the **last coefficient in `resultsNames(dds)`** -- which depends on factor level order and design formula order. Always specify the contrast explicitly. Tutorials that hard-code `results(dds)` are setting an example that breaks the moment another factor is added.

## Algorithmic Taxonomy

| Test / estimator | What it tests | When mandatory | Failure mode |
|------------------|---------------|----------------|--------------|
| Wald | One coefficient = 0 | Two-level factor or single contrast | Anti-conservative with many low-count outliers |
| LRT (`test='LRT'`, `reduced=`) | Joint effect of dropped terms (>=1 df) | Multi-level factor, omnibus, interaction with >1 df | Reports LFC of the LAST coefficient, not omnibus -- read p-value but never report the LFC as "the effect" |
| `lfcShrink(type='apeglm')` (Zhu 2019) | Posterior LFC under heavy-tailed Cauchy prior | DEFAULT for ranking and visualization | Requires `coef=`; cannot use `contrast=` or numeric vectors |
| `lfcShrink(type='ashr')` (Stephens 2017) | Posterior LFC under unimodal prior; reports `lfsr`/`svalue` with `svalue=TRUE` | Arbitrary contrasts via `contrast=` | Slightly different inferential frame (sign-error rather than null-FDR) |
| `lfcShrink(type='normal')` | Posterior LFC under zero-centered normal; accepts `coef=` or `contrast=` (with `res=`) | Quasi-deprecated since v1.16; only path to get shrunken p-values | Cannot be used with formulas containing interaction terms |
| TREAT / `lfcThreshold=` (McCarthy & Smyth 2009) | LFC magnitude exceeds threshold tau | Want FDR control for "|LFC| > 1.5x" claims | Conservative; use only when threshold is biologically pre-specified |

## Decision Tree by Scenario

| Scenario | Recommended path | Why |
|----------|------------------|-----|
| Two-group bulk RNA-seq, n>=3/group | `DESeq()` + `results(name=...)` + `lfcShrink(coef=..., type='apeglm')` | Modern default; apeglm is the right prior |
| Factor with 3+ levels, "any change" question | `DESeq(test='LRT', reduced=~1)`; read padj only | Wald + p-value combining is wrong |
| Interaction `~ genotype * treatment` | Build combined factor `group = paste(genotype, treatment)` and design `~ 0 + group`; contrast pairs of interest | Avoids the resultsNames trap; works with apeglm via relevel |
| Paired design (tumor/normal same patient) | `~ patient + tissue`; pairing variable FIRST | Absorbs subject variability; n_paired effective sample size |
| Salmon/kallisto input | `tximport()` -> `DESeqDataSetFromTximport()` | Carries length offsets automatically; `DESeqDataSetFromMatrix(round(...))` loses length correction |
| n=2/group, no choice | Continue but report results as exploratory; consider edgeR QL F-test as sensitivity | Schurch 2016 *RNA* 22:839: all tools miss 20-40% of true positives at n=3 |
| Single-cell pseudobulk (counts aggregated per donor) | DESeq2 standard pipeline on pseudobulk matrix | Crowell 2020 *Nat Commun* 11:6077: pseudobulk avoids the FDR inflation of cell-level DE |
| Many DE genes expected (>50% of genome) | `estimateSizeFactors(controlGenes=stable)` or spike-in normalization | Median-of-ratios assumes most genes unchanged |
| GSEA preranked input | Shrunken LFC OR `stat` (Wald Z) as the rank | Unshrunken LFC dominated by low-count noise |
| Cross-sample heatmap, PCA, ML feature | `vst(dds, blind=FALSE)` (or `rlog` if n<30 and library sizes vary >4x) | Raw counts make PC1 = library size |

## Standard Workflow

**Goal:** Take a raw integer count matrix and a sample table to a ranked, shrunken DE result table.

**Approach:** Construct DESeqDataSet with the design formula, set reference levels explicitly, run the pipeline, extract by explicit contrast, shrink for downstream use.

```r
library(DESeq2)
library(apeglm)

dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
dds$condition <- relevel(dds$condition, ref = 'control')

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds, name = 'condition_treated_vs_control', alpha = 0.05)
res_shrunk <- lfcShrink(dds, coef = 'condition_treated_vs_control', type = 'apeglm')

summary(res)
sig <- subset(res, padj < 0.05)
```

The reference level fix is non-cosmetic: DESeq2 picks alphabetically if not told otherwise, so `c('Treated','Untreated')` makes 'Treated' the reference and the LFC reads inverted. Set it BEFORE `DESeq()`.

## Tximport (Salmon / kallisto / RSEM)

For salmon/kallisto/RSEM input, use `DESeqDataSetFromTximport()` which carries the per-sample length matrix as an offset automatically:

```r
library(tximport)
txi <- tximport(files, type = 'salmon', tx2gene = tx2gene)
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ condition)
dds <- DESeq(dds)
```

The `tximport(..., countsFromAbundance='lengthScaledTPM')` form is for limma-voom (no offset mechanism). Full mechanics, the four `countsFromAbundance` options, RSEM zero-length traps, and tximeta provenance: see `expression-matrix/counts-ingest`.

## Design Formulas and the resultsNames Trap

**Goal:** Encode batch, paired, and interaction structure correctly and extract the intended contrast.

**Approach:** Put the variable of interest LAST for readability, but never trust the default `results(dds)` -- inspect `resultsNames(dds)` and pass `name=` or `contrast=` explicitly.

```r
design(dds) <- ~ batch + condition
dds <- DESeq(dds)
resultsNames(dds)
# "Intercept" "batch_B_vs_A" "condition_treated_vs_control"

res <- results(dds, name = 'condition_treated_vs_control')
```

Interaction design with the canonical trap:

```r
design(dds) <- ~ genotype + treatment + genotype:treatment
dds <- DESeq(dds)
resultsNames(dds)
# "Intercept" "genotype_KO_vs_WT" "treatment_drug_vs_vehicle" "genotypeKO.treatmentdrug"
```

- `results(dds, name='treatment_drug_vs_vehicle')` returns the drug effect IN THE WT REFERENCE only, NOT a marginal average. This is the single most common misinterpretation.
- `results(dds, name='genotypeKO.treatmentdrug')` returns the DIFFERENCE in drug effect between KO and WT.
- Drug effect in KO requires summing: `results(dds, contrast=list(c('treatment_drug_vs_vehicle','genotypeKO.treatmentdrug')))`.

Cleaner alternative for interactions when many pairwise contrasts are needed: combined factor + `~ 0 + group`.

```r
dds$group <- factor(paste(dds$genotype, dds$treatment, sep = '_'))
design(dds) <- ~ 0 + group
dds <- DESeq(dds)
res_drug_in_ko <- results(dds, contrast = c('group', 'KO_drug', 'KO_vehicle'))
```

## Wald vs LRT

**Goal:** Choose Wald for single-coefficient hypotheses, LRT for joint hypotheses involving more than 1 df.

**Approach:** Wald is default. LRT is mandatory for multi-level factors tested as "any change", interactions with >1 df, and ANOVA-style omnibus tests. With LRT, the reported LFC is for the LAST coefficient in `resultsNames(dds)` -- not the omnibus effect.

```r
dds <- DESeq(dds, test = 'LRT', reduced = ~ batch)
res_lrt <- results(dds)
# res_lrt$padj is the LRT joint p-value (correct).
# res_lrt$log2FoldChange is for the last coefficient in resultsNames -- NOT an omnibus summary.
```

When reporting LRT results, name what the LFC actually represents or extract specific Wald coefficients per level for the effect-size table.

## LFC Shrinkage -- Three Flavors, Three Failure Modes

**Goal:** Get a stable effect-size estimate appropriate for ranking, GSEA input, volcano-plot x-axis, and reporting.

**Approach:** Default to apeglm. Switch to ashr when arbitrary contrasts are needed. Never use normal for new analyses.

```r
res_apeglm <- lfcShrink(dds, coef = 'condition_treated_vs_control', type = 'apeglm')
res_ashr   <- lfcShrink(dds, contrast = c('condition','treated','control'), type = 'ashr')
```

| Method | Prior | Accepts | Use when |
|--------|-------|---------|----------|
| apeglm | Cauchy (heavy-tailed) | `coef=` only | Default; preserves large effects, suppresses low-count noise |
| ashr | Unimodal scale-mixture | `coef=` or `contrast=` (incl. numeric) | Need contrast= for interaction sums or pairwise from `~ 0 + group` |
| normal | Zero-centered normal | `coef=` only; not interaction designs | Only when the old shrunken p-value is required (legacy) |

The apeglm-cannot-use-contrast footgun: if the question is "drug effect in KO" from `~ genotype * treatment`, apeglm cannot directly shrink that contrast. Workarounds: (a) rebuild as combined factor `~ 0 + group` and relevel so the desired comparison is a coefficient; (b) use ashr; (c) accept the unshrunken LFC for that one comparison.

p-values do NOT change when shrinking. `lfcShrink()` preserves the Wald p-value from `results()`. See the Single Most Important Insight at the top.

## Independent Filtering, Cook's Outliers, padj=NA

`padj = NA` has three distinct causes (independent filtering, Cook's outlier, all-zero in a group), each with a different remediation -- see `de-results` for the full diagnostic table, IHW alternative, and recovery code.

Two DESeq2-specific points worth knowing at this layer:

- Cook's distance filtering is NOT computed for continuous covariates -- a continuous-covariate analysis has effectively no automatic outlier filtering. Disable Cook's only when the outlier IS the signal (`results(dds, cooksCutoff = FALSE)`). At n>=7 per group, `DESeq()` REPLACES outliers via `replaceOutliers()` and refits (`minReplicatesForReplace = 7`).
- Pre-filtering (`rowSums(counts(dds)) >= 10`) is for memory and speed ONLY. It does NOT replace independent filtering, which operates downstream at `results()` time. Independent filtering can be swapped for IHW via `results(dds, filterFun = ihw)`.

## VST vs rlog vs normTransform (Visualization Only)

**Goal:** Produce homoskedastic log2-scale counts for PCA, heatmaps, clustering, ML features.

**Approach:** Use `vst()` by default. Switch to `rlog()` only for n<30 with size factors varying >4x. Never use `normTransform()` for distance-based plots. Never use VST/rlog values as input to DE.

```r
vsd <- vst(dds, blind = FALSE)
rld <- rlog(dds, blind = FALSE)
```

The `vst()` function default is `blind=TRUE`, but the current DESeq2 vignette recommends `blind=FALSE` for any downstream visualization AFTER the model is fit (it uses the design when fitting dispersions; appropriate when the design is already settled). Reserve `blind=TRUE` for unsupervised QC where the design should not influence the transformation (e.g., "is this sample consistent with its group?"). The vignette has flip-flopped over the years on which to recommend by default -- pass `blind=` explicitly.

`vst()` uses 1000 most-variable genes to fit the dispersion trend by default. With <1000 genes after filtering, set `nsub` lower.

## betaPrior Deprecation Timeline

DESeq2 v1.0-v1.15: `betaPrior = TRUE` was the default. Shrinkage was baked in and p-values were computed on the shrunken estimates. v1.16 (2017) flipped the default to `FALSE` and introduced `lfcShrink()`. Today: `betaPrior = TRUE` is quasi-deprecated and is the ONLY path to a p-value of the shrunken estimate. Most users do not want this. `lfcShrink()` with apeglm and the unshrunken Wald p-value is the modern compromise.

The contrast= vs name= numerical difference that older tutorials describe was specific to `betaPrior = TRUE`. With the current default, `name=` and `contrast=` for the same comparison return identical LFCs.

## Size Factor Alternatives

`estimateSizeFactors()` defaults to `type='ratio'` (median-of-ratios). Edge cases:

| Situation | Use |
|-----------|-----|
| Zero counts in some samples for many genes | `type='poscounts'` -- uses only positive entries per gene |
| Very small libraries, hard to converge | `type='iterate'` |
| Spike-ins (ERCC) or known stable housekeeping genes | `controlGenes = indices` |
| Majority-DE biology (prokaryotic stress, viral host shutoff, MYC amplification) | `controlGenes` with curated stable genes; or spike-in SBN (Jiang 2011 *Genome Res* 21:1543) |

Single-cell pseudobulk: most genes have zeros across donors, so `type='poscounts'` is often required for pseudobulk DESeq2.

## Per-Method Failure Modes

### apeglm refuses arbitrary contrasts

**Trigger:** Question of the form "drug effect in KO genotype" from `~ genotype * treatment`; user tries `lfcShrink(dds, contrast=list(...), type='apeglm')`.

**Mechanism:** apeglm fits per-coefficient priors. A numeric or list contrast is a linear combination of coefficients, not a coefficient -- no prior to apply.

**Symptom:** Error: "type='apeglm' shrinkage only for use with 'coef'"

**Fix:** Rebuild design as `~ 0 + group` with `group = paste(genotype, treatment)`, relevel so the comparison is a coefficient, refit. Or use `type='ashr'` which accepts contrasts.

### LRT reports the wrong LFC

**Trigger:** Multi-level factor analyzed with `DESeq(dds, test='LRT', reduced=~1)`; user reports the `log2FoldChange` column as "the effect".

**Mechanism:** The LRT p-value is the omnibus test of "any difference among levels". The LFC reported by `results()` after LRT is for the LAST coefficient in `resultsNames(dds)`, which is one specific level-vs-reference comparison.

**Symptom:** A 4-level factor produces a single LFC value per gene; reviewer asks "the effect of which condition?"

**Fix:** Treat LRT padj as a screen for "any change". For effect sizes, extract per-level Wald coefficients individually via `results(dds, name='<specific coefficient>')` for each non-reference level.

### Cook's silences the gene of interest

**Trigger:** Rare-disease cohort or CNV-amplified patient where ONE sample drives the biology; that gene has `padj=NA`.

**Mechanism:** Cook's distance flagged the patient as an outlier and zeroed the gene's p-value (or replaced its counts if n>=7).

**Symptom:** A gene of clear biological interest comes back as NA in the results table even though the count matrix shows the expected pattern.

**Fix:** `results(dds, cooksCutoff = FALSE)`. Optionally cross-validate the result by running a sensitivity analysis with and without the outlier sample.

### Independent filtering kills the master regulator

**Trigger:** A transcription factor expressed at ~10 counts but consistently across all samples shows `padj=NA` despite obvious biology.

**Mechanism:** baseMean is below the data-driven independent filtering threshold; the gene was excluded from FDR adjustment.

**Symptom:** Low-count genes with clean signal end up NA.

**Fix:** `results(dds, independentFiltering = FALSE)`; or `filterFun = ihw` from the IHW package (often less aggressive on low-count genes).

### Parametric dispersion-mean trend doesn't fit the cloud

**Trigger:** `plotDispEsts(dds)` shows the red parametric trend curve nowhere near the cloud of gene-wise (blue) and final (black) dispersion estimates.

**Mechanism:** Default `fitType='parametric'` assumes `dispersion ~ a/mean + b`. Fails when the experiment has very few samples per group with highly heterogeneous biology, many very-low-count genes pulling the trend, or a continuous covariate driving large variability.

**Symptom:** Curved trend that doesn't match the cloud; mismatch shows up most clearly in low-baseMean genes.

**Fix:** Refit with `DESeq(dds, fitType='local')` (local regression) or `fitType='mean'` (flat trend); compare `plotDispEsts` between fits and pick the one tracking the cloud. Falling back to `'mean'` is a sign the data is unusual; investigate before trusting results.

### Median-of-ratios fails on prokaryotic stress

**Trigger:** Bacterial RNA-seq under stress where >50% of genes change in one direction; PCA shows huge global shift.

**Mechanism:** Median-of-ratios assumes most genes are not DE. Under massive global perturbation, the reference is dominated by DE genes; size factors absorb the biology.

**Symptom:** MA plot shows the bulk cloud shifted off zero; reported fold changes don't match qPCR.

**Fix:** Use `controlGenes` with curated stable housekeeping genes; or spike-in normalization (Jiang 2011); or RUVg with negative controls.

## PyDESeq2 (Python alternative)

```python
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

dds = DeseqDataSet(counts=count_df, metadata=metadata, design='~condition')
dds.deseq2()
stat_res = DeseqStats(dds, contrast=('condition', 'treated', 'control'))
stat_res.summary()
results_df = stat_res.results_df
```

PyDESeq2 0.5+ supports Wald, multi-factor designs, and apeglm shrinkage. No LRT yet. Results are numerically close to R DESeq2 but with small differences from MLE solver choice.

## Prokaryotic RNA-seq

- Non-spliced aligners (BWA-MEM, Bowtie2) -- no introns.
- Polycistronic operons cause read-through between adjacent genes; confirm gene boundaries in the GFF.
- rRNA depletion essential (80-95% rRNA without poly-A selection, which prokaryotes lack anyway).
- Median-of-ratios fails under stress (see failure mode above) -- use `controlGenes` or spike-ins.
- KEGG organism codes are strain-specific (e.g., `pae` for P. aeruginosa PAO1): `clusterProfiler::search_kegg_organism()`.
- Annotation comes from Prokka or Bakta GFF; Ensembl/biomaRt are eukaryote-only.

## Common errors

| Error / symptom | Cause | Fix |
|-----------------|-------|-----|
| `design matrix not full rank` | Confounded covariates | `alias(model.matrix(design, coldata))$Complete` to find the redundant column |
| `counts matrix should be integers` | Salmon/kallisto/RSEM counts are fractional | Use `DESeqDataSetFromTximport()` -- it rounds AND carries length offsets |
| Wrong sign of LFC vs expected | Reference level set alphabetically | `relevel()` BEFORE `DESeq()` |
| `padj = NA` for biologically meaningful gene | Independent filtering or Cook's outlier | See padj=NA section |
| LRT LFC doesn't match Wald LFC for the same comparison | LRT reports last coefficient; Wald reports the named coefficient | Extract specific Wald per level for the effect size |
| `summary(res)` shows fewer DE genes than expected | `summary()` default `alpha=0.1`, NOT the alpha passed to `results()` | `summary(res, alpha = 0.05)` |

## References

- Anders S, Huber W. 2010. Differential expression analysis for sequence count data. *Genome Biol* 11(10):R106. doi:10.1186/gb-2010-11-10-r106
- Love MI, Huber W, Anders S. 2014. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol* 15(12):550. doi:10.1186/s13059-014-0550-8
- Zhu A, Ibrahim JG, Love MI. 2019. Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences. *Bioinformatics* 35(12):2084-2092. doi:10.1093/bioinformatics/bty895
- Stephens M. 2017. False discovery rates: a new deal. *Biostatistics* 18(2):275-294. doi:10.1093/biostatistics/kxw041
- Ignatiadis N, Klaus B, Zaugg JB, Huber W. 2016. Data-driven hypothesis weighting increases detection power in genome-scale multiple testing. *Nat Methods* 13(7):577-580. doi:10.1038/nmeth.3885
- Bourgon R, Gentleman R, Huber W. 2010. Independent filtering increases detection power for high-throughput experiments. *PNAS* 107(21):9546-9551. doi:10.1073/pnas.0914005107
- McCarthy DJ, Smyth GK. 2009. Testing significance relative to a fold-change threshold is a TREAT. *Bioinformatics* 25(6):765-771. doi:10.1093/bioinformatics/btp053
- Soneson C, Love MI, Robinson MD. 2015. Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences. *F1000Res* 4:1521. doi:10.12688/f1000research.7563.2
- Schurch NJ et al. 2016. How many biological replicates are needed in an RNA-seq experiment and which differential expression tool should you use? *RNA* 22(6):839-851. doi:10.1261/rna.053959.115
- Crowell HL et al. 2020. muscat detects subpopulation-specific state transitions from multi-sample multi-condition single-cell transcriptomics data. *Nat Commun* 11:6077. doi:10.1038/s41467-020-19894-4
- Jiang L et al. 2011. Synthetic spike-in standards for RNA-seq experiments. *Genome Res* 21(9):1543-1551. doi:10.1101/gr.121095.111

## Related Skills

- edger-basics - Cross-check or use when n<5/group; QL F-test framework
- de-results - padj=NA handling, IHW, TREAT, GSEA input preparation, gene annotation
- de-visualization - MA, volcano (with shrunken LFC), PCA, heatmap, dispersion plot
- batch-correction - Include batch in design vs Nygaard 2016 cardinal sin
- timeseries-de - DESeq2 LRT with splines for time-course
- expression-matrix/counts-ingest - tximport, featureCounts, STAR output decisions
- expression-matrix/normalization - RLE/TMM/VST/rlog mechanics and failure modes
- expression-matrix/metadata-joins - Reference level, paired design, interaction parameterization
- expression-matrix/gene-id-mapping - Annotating DE results with symbols
- rna-quantification/tximport-workflow - Detailed tximport mechanics
- pathway-analysis/gsea - Ranked-list input from DE
- pathway-analysis/go-enrichment - ORA with proper background
- data-visualization/volcano-and-ma-plots - Custom volcano with apeglm-shrunken LFC

---
name: bio-microbiome-differential-abundance
description: Tests which individual taxa differ between groups on an amplicon ASV/feature table (phyloseq) using compositionally-aware methods - ALDEx2 (Dirichlet-MC CLR, conservative), ANCOM-BC2/ANCOMBC (sampling-fraction bias correction, structural zeros, passed_ss, default p_adj_method=holm), MaAsLin2/MaAsLin3 (multivariable GLM, random effects, prevalence/abundance split), LinDA (CLR mixed-model regression), ZicoSeq (permutation FDR), LEfSe, and q2-composition ancombc. Covers why the hit list depends more on the DA tool than the biology (Nearing benchmark) so the deliverable is a CONSENSUS of >=2 tools, why a relative change is not absolute without a load anchor, the prevalence-filter knob, BH/FDR plus an effect-size floor, and why DESeq2/edgeR misfire here. Use when finding differentially abundant taxa, handling covariates or longitudinal designs, or choosing a method. Whole-community diversity -> diversity-analysis; shotgun DA -> metagenomics/metagenome-visualization; CoDA theory -> metagenomics/abundance-estimation
tool_type: r
primary_tool: ALDEx2
---

## Version Compatibility

Reference examples tested with: ALDEx2 1.34+, ANCOMBC 2.4+, Maaslin2 1.16+, MicrobiomeStat 1.2+ (LinDA), GUniFrac 1.8+ (ZicoSeq), phyloseq 1.46+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

ANCOM-BC2 changed argument names between `ancombc()` and `ancombc2()`, and its default `p_adj_method` is `holm`, not `BH` - confirm both against the installed version. The MaAsLin3 `maaslin3()` API differs from MaAsLin2's `Maaslin2()`.

# Differential Abundance Testing

**"Find which taxa differ between my groups"** -> Run two or more compositionally-aware DA tools and report their consensus - because the significant-taxa list is a property of the tool as much as of the sample, and a relative-abundance change is not an absolute change.
- R: `ALDEx2::aldex(counts, conds, test='t', effect=TRUE, denom='all')` then a second tool (`ANCOMBC::ancombc2()` or `MicrobiomeStat::linda()`)

Scope: per-taxon DA on an amplicon feature table. Whole-community alpha/beta/PERMANOVA -> diversity-analysis. Shotgun profiler-table DA -> metagenomics/metagenome-visualization. Shared compositional/closure/CLR/zero theory -> metagenomics/abundance-estimation. Collapse ASVs to genus/species first -> taxonomy-assignment. QIIME2 CLI route -> qiime2-workflow.

## The Single Most Important Modern Insight -- Which Taxa Are "Significant" Depends More on the Tool Than on the Biology

Run ALDEx2, ANCOM-BC2, MaAsLin2, and LinDA on the same ASV table and the four significant-taxa lists overlap but disagree (Nearing 2022 *Nat Commun* 13:342, across 38 datasets). So the deliverable is NOT "the differential taxa" - it is the CONSENSUS of >=2 compositionally-aware tools, every tool NAMED: the intersection is high-confidence, the union is exploratory, and a single-tool hit is tentative. Picking the tool with the prettiest volcano is p-hacking by software (uncorrected multiplicity hidden in the method menu). Three corollaries:

1. **A relative-abundance increase is not an absolute increase.** Microbiome counts are compositional - the sequencer fixes the total, so one taxon blooming forces every other taxon's proportion down (the blooming-taxon illusion). "Taxon X increased" is a statement about its SHARE unless an external load anchor (spike-in / flow cytometry / qPCR, see metagenomics/abundance-estimation) or MaAsLin3's absolute-abundance mode licenses an absolute claim.
2. **Uncorrected Wilcoxon/t-test on raw relative abundances is wrong twice in one line** - closure (reference-frame) AND multiple testing. But a BH-corrected simple test, honestly labelled as relative, can replicate BETTER than a sophisticated model (Pelto 2025): the forbidden thing is the uncorrected, closure-blind form, not simple tests per se.
3. **There is no settled best tool.** The benchmarks optimize different criteria, so they rank tools differently. Consensus-of-tools is the only stance that survives all of them.

## The Benchmark Landscape (no settled winner)

| Benchmark | Optimized for | Verdict |
|-----------|---------------|---------|
| Nearing 2022 *Nat Commun* 13:342 | cross-method consistency | ALDEx2 + ANCOM-II most consistent and most conservative; LEfSe/edgeR flag far more, agree less |
| Yang & Chen 2022 *Microbiome* 10:130 | FDR-power balance | ZicoSeq / LinDA / ANCOM-BC-family best |
| Yang & Chen 2023 *Brief Bioinform* 24:bbac607 | correlated (repeated-measures) designs | use a mixed-model-capable tool (LinDA, MaAsLin2, ANCOM-BC2) |
| Pelto 2025 *Brief Bioinform* 26(2):bbaf130 | cross-study replicability | elementary BH-corrected methods most replicable; ANCOM-BC2 worst |

Report the disagreement AS the result; verify current best practice against the latest tool docs rather than hard-coding one method.

## Tool Taxonomy

| Tool | Citation | Mechanism / role | When |
|------|----------|------------------|------|
| ALDEx2 | Fernandes 2014 *Microbiome* 2:15 | Dirichlet Monte-Carlo posterior + CLR; tests each draw; reports expected effect + BH-adjusted p | conservative two-group anchor; small-to-moderate n |
| ANCOM-BC2 | Lin & Peddada 2024 *Nat Methods* 21:83 | estimates per-sample sampling fraction and bias-corrects; structural zeros; pseudo-count sensitivity (`passed_ss`) | interpretable LFC + CI; covariates; multi-group |
| MaAsLin2 | Mallick 2021 *PLoS Comput Biol* 17:e1009442 | general (mixed) linear model on transformed abundance | multivariable / longitudinal / metadata-rich |
| MaAsLin3 | Nickols 2026 *Nat Methods* 23:554 | splits abundance (level when present) from prevalence (present/absent); absolute-abundance mode | prevalence-vs-abundance separation; load data available |
| LinDA | Zhou 2022 *Genome Biol* 23:95 | CLR regression with mode-based bias correction; asymptotic FDR | large cohorts; fast; native mixed model |
| ZicoSeq | Yang & Chen 2022 *Microbiome* 10:130 | reference-taxa normalization + permutation FDR; winsorization | covariates; non-parametric permutation p; strong FP control |
| LEfSe | Segata 2011 *Genome Biol* 12:R60 | Kruskal-Wallis + LDA effect size | exploratory biomarker discovery; NOT a formal FDR-controlled test |
| DESeq2 | Love 2014 *Genome Biol* 15:550 | RNA-seq median-of-ratios size factor | caveat only; geometric-mean reference dies on sparse zero-heavy tables |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Two groups, want a trustworthy conservative anchor | ALDEx2 | Dirichlet-MC + CLR; most reproducible/conservative (Nearing); gate on effect size |
| Need interpretable LFC + CI, structural zeros, multi-group | ANCOM-BC2 | models and corrects per-sample sampling fraction; global/pairwise/Dunnett/trend; `passed_ss` |
| Large cohort, covariates, speed, mixed model | LinDA | CLR regression + bias mode; asymptotic FDR; fast; random effects in the formula |
| Covariates + permutation-grounded non-parametric p | ZicoSeq | reference-taxa frame + permutation FDR |
| Longitudinal / many covariates / flexible GLM | MaAsLin2 | `fixed_effects` + `random_effects`; normalization/transform menu |
| Prevalence-vs-abundance separation or absolute abundance | MaAsLin3 | logistic prevalence model + abundance model; load-data hook |
| Inside a QIIME2 CLI pipeline | qiime composition ancombc + tabulate/da-barplot | native artifact flow (v1 ANCOM-BC; go to R for v2 `passed_ss`/multi-group) -> qiime2-workflow |
| Repeated / paired samples | any tool above WITH a random effect | ignoring subject structure is pseudo-replication |
| ALWAYS | run >=2 of the above, report the consensus | tool choice drives the hit list more than biology (Nearing 2022) |
| Shotgun species table, not amplicon | -> metagenomics/metagenome-visualization | same CoDA theory; different upstream pipeline |
| Uncorrected t-test/Wilcoxon on TSS proportions | DO NOT | closure biases the test and there is no FDR control |

## Filter Before Testing (a modeling knob, not housekeeping)

**Goal:** Drop rare features before testing so the BH denominator is not crushed and log/CLR transforms are well-behaved.

**Approach:** Keep features present in at least 10-25% of samples (and optionally a mean-abundance floor); declare the threshold and confirm the headline result is not knife-edge-sensitive to it. Every tool exposes this (`prv_cut`, `min_prevalence`, `prev.filter`).

```r
library(phyloseq)
ps <- readRDS('phyloseq_object.rds')
# prv_cut 0.10: a feature must appear in >= 10% of samples; raising to 0.25 removes more tests
# (smaller BH correction, more power on survivors) but discards rare-but-real taxa - a declared choice
keep <- filter_taxa(ps, function(x) sum(x > 0) >= 0.10 * nsamples(ps), TRUE)
```

## ALDEx2: The Conservative Floor of the Consensus

**Goal:** Identify taxa that differ between two groups while propagating the sampling uncertainty of low-count features.

**Approach:** Draw `mc.samples` Monte-Carlo instances from a Dirichlet posterior of the counts (this IS the zero handling - no explicit pseudocount), CLR-transform each instance against the geometric mean of all features (`denom='all'`), run the test on every draw, and report the EXPECTED effect size and BH-adjusted p over the draws.

```r
library(ALDEx2)
counts <- as.matrix(otu_table(ps))            # integer counts, taxa in ROWS
if (!taxa_are_rows(ps)) counts <- t(counts)
groups <- as.character(sample_data(ps)$Group)

# mc.samples 128: standard Monte-Carlo draws; 256+ for publication (more stable expected p)
res <- aldex(counts, groups, mc.samples = 128, test = 't', effect = TRUE, denom = 'all')
# we.eBH = Welch expected BH-adjusted p (report this, NOT we.ep); wi.eBH = Wilcoxon equivalent
# effect = median standardized effect = median(diff.btw / max(diff.win)); the primary decision variable
hits <- res[res$we.eBH < 0.05 & abs(res$effect) > 1, ]   # q AND effect floor (Gloor: gate on effect, not p alone)
```

Gate on effect size AND q, not p alone: with large n trivially small CLR differences become "significant," and Gloor's own guidance is that `|effect| > 1` is a strong ~2-SD signal. For >2 groups use `aldex.kw()`; for covariates the `aldex.glm()` + model.matrix route works but ALDEx2 is weakest here - prefer ANCOM-BC2/LinDA/MaAsLin2 for serious covariate or random-effect modeling.

## ANCOM-BC2: Bias-Corrected LFC With a Sensitivity Safeguard

**Goal:** Estimate an interpretable bias-corrected log-fold-change per taxon, with covariate adjustment, structural-zero handling, and a flag for hits that are hostage to the pseudo-count.

**Approach:** Model log(observed count) as a function of covariates, estimate each sample's log sampling fraction as an offset and subtract it, then refit across a range of pseudo-counts and record how often each q-value flips (`passed_ss`).

```r
library(ANCOMBC)
out <- ancombc2(data = ps, fix_formula = 'Group + Age + Sex',
                rand_formula = NULL,        # '(1 | SubjectID)' for repeated measures - see Failure Modes
                p_adj_method = 'BH',        # DEFAULT is 'holm'; set 'BH' deliberately for FDR
                prv_cut = 0.10, lib_cut = 1000,
                group = 'Group', struc_zero = TRUE, pseudo_sens = TRUE,
                global = FALSE, pairwise = FALSE, n_cl = 2)
res <- out$res
# a confident hit is BOTH significant AND robust to the pseudo-count. ANCOM-BC2 suffixes the
# diff_/passed_ss_ columns with the literal model-matrix coefficient (variable + factor level,
# verbatim case, e.g. 'Grouptreated') - match it by pattern rather than hard-coding the case.
dcol <- grep('^diff_Group', names(res), value = TRUE)[1]
robust <- res[res[[dcol]] & res[[sub('^diff_', 'passed_ss_', dcol)]], ]
```

`passed_ss` is the most valuable ANCOM-BC2-specific feature: a CLR/log model on sparse data is hostage to the zero-replacement constant, and `passed_ss` quantifies that per taxon. A hit with `passed_ss == FALSE` depends on the arbitrary pseudo-count - do not report it as confident. For >2 groups set `global=TRUE` (omnibus), `pairwise=TRUE` (mdFDR-controlled pairs), `dunnet=TRUE`, or `trend=TRUE`; results land in `out$res_global`/`res_pair`/`res_dunn`/`res_trend`.

## LinDA: Fast CLR Regression With Native Mixed Models

**Goal:** Get FDR-controlled log2-fold-changes on a large cohort, including repeated-measures designs, without Monte-Carlo or EM cost.

**Approach:** Fit ordinary linear regression on the CLR-transformed table covariate by covariate, estimate the compositional bias as the mode of the per-feature coefficients and subtract it; a random effect in the formula makes it a linear mixed model.

```r
library(MicrobiomeStat)
otu <- as.data.frame(otu_table(ps)); if (!taxa_are_rows(ps)) otu <- t(otu)
meta <- as.data.frame(sample_data(ps))
fit <- linda(feature.dat = otu, meta.dat = meta,
             formula = '~ Group + Age + (1 | SubjectID)',   # random effect -> mixed model
             feature.dat.type = 'count', prev.filter = 0.10, alpha = 0.05)
fit$output[[1]]   # names(fit$output) are the model-matrix coefficient columns (e.g. 'Grouptreated' - the factor level keeps its case); per-feature: log2FoldChange, lfcSE, stat, pvalue, padj, reject
```

LinDA is the natural fast modern entry in a consensus panel and the cleanest route to mixed models. Yang & Chen rate it among the best FDR-power trade-offs.

## MaAsLin2 / MaAsLin3 and ZicoSeq (the rest of the panel)

**Goal:** Fit covariate-rich or longitudinal differential-abundance models, or add a permutation-based panel member, when ALDEx2/ANCOM-BC2/LinDA do not cover the design.

**Approach:** Use MaAsLin2/3 for multivariable GLMs with random effects, or ZicoSeq for a non-parametric permutation-FDR test against empirically selected reference taxa.

MaAsLin2 fits a flexible per-feature GLM; its package DEFAULT is TSS + LOG + LM (not CLR), and `random_effects` is the canonical route for longitudinal designs. NOTE the orientation gotcha: it expects features in COLUMNS, samples in rows.

```r
library(Maaslin2)
fit <- Maaslin2(input_data = as.data.frame(t(otu)), input_metadata = meta,
                output = 'maaslin2_out', fixed_effects = c('Group', 'Age'),
                random_effects = c('SubjectID'),
                normalization = 'TSS', transform = 'LOG', analysis_method = 'LM',
                min_prevalence = 0.10, max_significance = 0.05)
# writes all_results.tsv / significant_results.tsv with columns feature, metadata, coef, pval, qval
```

MaAsLin3 (`maaslin3()`) splits each feature into an abundance model (level when present) and a logistic prevalence model (present/absent) tested jointly, and can ingest total-load measurements for absolute-abundance inference. ZicoSeq (`GUniFrac::ZicoSeq()`) winsorizes, posterior-samples, normalizes against empirically selected reference taxa, and returns permutation FDR (`zc$p.adj.fdr`) - a non-parametric panel member that accepts covariates via `adj.name`.

## Consensus: Intersect the Tools

**Goal:** Convert two or more per-tool hit sets into a confidence-graded result instead of one tool's answer.

**Approach:** Collect the significant feature SETS (BH within each tool), then report the intersection as high-confidence, the union as exploratory, and tabulate, per taxon, how many of N tools agree and which ones. Never pool p-values across tools.

```r
sig_aldex <- rownames(res)[res$we.eBH < 0.05 & abs(res$effect) > 1]
sig_linda <- rownames(fit$output[[1]])[fit$output[[1]]$reject]   # [[1]] = the group coefficient (named 'Grouptreated')
confident  <- intersect(sig_aldex, sig_linda)   # high-confidence
exploratory <- union(sig_aldex, sig_linda)       # report with the tool that found each
```

## Per-Method Failure Modes

### Cherry-picking the tool with the prettiest result
**Trigger:** running several tools and reporting only the one(s) that flag the favored taxon. **Mechanism:** that is uncorrected multiplicity hidden in the method menu (p-hacking by software). **Symptom:** "the recommended method found X" with no mention of the tools that disagreed. **Fix:** decide the panel a priori, report ALL tools, intersect for confident hits, disclose disagreement.

### Uncorrected Wilcoxon/t-test on relative abundances
**Trigger:** a per-taxon Wilcoxon/t-test on TSS proportions with no FDR correction. **Mechanism:** closure makes a naive test call every taxon "decreased" when one blooms, and hundreds of uncorrected tests inflate false positives. **Symptom:** dozens of "significant" taxa, all in the same direction, no q-values. **Fix:** use a CoDA/reference-frame tool; if a simple test is used, BH-correct it and label the comparison as relative (Pelto 2025).

### Pseudo-replication of repeated measures
**Trigger:** longitudinal/paired samples treated as independent rows. **Mechanism:** fewer independent units than rows inflates significance. **Symptom:** implausibly small p-values on a small subject count. **Fix:** a random effect - ANCOM-BC2 `rand_formula='(1|SubjectID)'`, MaAsLin2 `random_effects='SubjectID'`, LinDA `(1|SubjectID)` in the formula. Cross-check ANCOM-BC2 mixed-model output against LinDA/MaAsLin2 (GitHub issue #111 reported `rand_formula` correctness problems in some versions).

### Prevalence filter set blindly
**Trigger:** an undeclared prevalence cut, or none at all. **Mechanism:** the cut decides which taxa are even tested and thus the BH landscape - it is a modeling choice. **Symptom:** the hit list changes materially between `prv_cut=0.1` and `0.25`. **Fix:** declare and justify the threshold; confirm the headline result survives moving it.

### Relative change reported as absolute
**Trigger:** "taxon X doubled" from a closed table with no load data. **Mechanism:** one taxon blooming compresses every other proportion. **Symptom:** whole-community "depletion" that is really one taxon rising. **Fix:** anchor to load (spike-in/flow/qPCR) or MaAsLin3 absolute mode; otherwise state the claim is relative.

### DESeq2/edgeR on a sparse 16S table
**Trigger:** RNA-seq median-of-ratios / TMM on a zero-heavy ASV table. **Mechanism:** the geometric-mean size-factor reference collapses on zeros and the "most features unchanged" assumption is violated. **Symptom:** degenerate size factors, errors, or inflated hit counts that disagree with CoDA tools (Nearing). **Fix:** use a compositional tool; if DESeq2 is unavoidable, the `poscounts` estimator is the minimum mitigation - present as a caveat, not a recipe.

### ANCOM-BC2 hit held hostage by the pseudo-count
**Trigger:** reporting `diff_* == TRUE` without checking `passed_ss_*`. **Mechanism:** significance depends on the arbitrary zero-replacement constant. **Symptom:** a hit that vanishes when the pseudo-count changes. **Fix:** require `diff_* & passed_ss_*` for a confident call.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Prevalence cut 10-25% (`prv_cut`/`min_prevalence`/`prev.filter`) | Nearing 2022; tool defaults (0.10) | rare features carry little information and crush the BH denominator; declare the value and test sensitivity |
| BH q <= 0.05 across taxa, within each tool | Benjamini-Hochberg 1995 *JRSS B* 57:289 | hundreds-thousands of features make uncorrected p meaningless; do not pool p across tools |
| ALDEx2 `|effect| > 1` (with q <= 0.05) | Gloor 2016 *J Comput Graph Stat* 25:971 | effect is a standardized median-ratio; ~2 SD is a strong signal; large n makes trivial diffs "significant" |
| ALDEx2 `mc.samples` = 128 (256+ for publication) | Fernandes 2014 *Microbiome* 2:15 | Monte-Carlo draws; more draws stabilize the expected p |
| ANCOM-BC2 `passed_ss == TRUE` required | Lin & Peddada 2024 *Nat Methods* 21:83 | flags hits whose significance is hostage to the pseudo-count |
| Consensus of >=2 compositionally-aware tools | Nearing 2022 *Nat Commun* 13:342 | tool choice drives the hit list more than biology; intersection = confident |
| ZicoSeq permutations `perm.no` >= 99 | Yang & Chen 2022 *Microbiome* 10:130 | permutation FDR resolution; raise for finer tail p |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| ALDEx2 returns NA effects / errors | proportions or non-integer matrix passed | feed integer COUNTS with taxa in rows |
| `passed_ss` column missing | `pseudo_sens = FALSE` | set `pseudo_sens = TRUE` (the default) |
| Far fewer hits than expected | ANCOM-BC2 `p_adj_method` left at `holm` | set `p_adj_method = 'BH'` deliberately if FDR is wanted |
| MaAsLin2 finds nothing / orientation error | features in rows, not columns | transpose so samples are rows, features columns |
| Mixed-model hits disagree across tools | `rand_formula` correctness varies by version | cross-check ANCOM-BC2 against LinDA/MaAsLin2 |
| Tools disagree on the hit list | normal - tool choice drives results | report the consensus and the disagreement, do not cherry-pick |
| Many "depleted" taxa in a host/plant sample | host mitochondria/chloroplast 16S inflates the table | filter Mitochondria/Chloroplast features (see taxonomy-assignment) before DA |
| Contaminant ASVs among the hits (low-biomass) | reagent kitome not removed before DA | run decontam upstream with negative controls (amplicon-processing; metagenomics/contamination-controls) |

## References

- Fernandes AD, Reid JNS, Macklaim JM, McMurrough TA, Edgell DR, Gloor GB. 2014. Unifying the analysis of high-throughput sequencing datasets: characterizing RNA-seq, 16S rRNA gene sequencing and selective growth experiments by compositional data analysis. *Microbiome* 2:15.
- Gloor GB, Macklaim JM, Fernandes AD. 2016. Displaying variation in large datasets: plotting a visual summary of effect sizes. *J Comput Graph Stat* 25:971-979.
- Lin H, Peddada SD. 2020. Analysis of compositions of microbiomes with bias correction (ANCOM-BC). *Nat Commun* 11:3514.
- Lin H, Peddada SD. 2024. Multigroup analysis of compositions of microbiomes with covariate adjustments and repeated measures (ANCOM-BC2). *Nat Methods* 21:83-91.
- Mallick H, Rahnavard A, McIver LJ, et al. 2021. Multivariable association discovery in population-scale meta-omics studies. *PLoS Comput Biol* 17:e1009442.
- Nickols WA, Kuntz T, Shen J, et al. 2026. MaAsLin 3: refining and extending generalized multivariable linear models for meta-omic association discovery. *Nat Methods* 23:554-564.
- Zhou H, He K, Chen J, Zhang X. 2022. LinDA: linear models for differential abundance analysis of microbiome compositional data. *Genome Biol* 23:95.
- Yang L, Chen J. 2022. A comprehensive evaluation of microbial differential abundance analysis methods: current status and potential solutions. *Microbiome* 10:130.
- Yang L, Chen J. 2023. Benchmarking differential abundance analysis methods for correlated microbiome sequencing data. *Brief Bioinform* 24:bbac607.
- Pelto J, Auranen K, Kujala JV, Lahti L. 2025. Elementary methods provide more replicable results in microbial differential abundance analysis. *Brief Bioinform* 26(2):bbaf130.
- Nearing JT, Douglas GM, Hayes MG, et al. 2022. Microbiome differential abundance methods produce different results across 38 datasets. *Nat Commun* 13:342.
- Segata N, Izard J, Waldron L, et al. 2011. Metagenomic biomarker discovery and explanation. *Genome Biol* 12:R60.
- Love MI, Huber W, Anders S. 2014. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biol* 15:550.
- Benjamini Y, Hochberg Y. 1995. Controlling the false discovery rate: a practical and powerful approach to multiple testing. *J R Stat Soc Series B* 57:289-300.

## Related Skills

- diversity-analysis - Whole-community alpha/beta/PERMANOVA; answer "do the communities differ" before "which taxa differ"
- taxonomy-assignment - Collapse ASVs to genus/species before per-taxon testing
- amplicon-processing - Produces the ASV feature table tested here
- qiime2-workflow - The qiime composition ancombc CLI route
- metagenomics/abundance-estimation - Shared compositional/closure/CLR/zero/load-anchor theory
- metagenomics/metagenome-visualization - The same DA mechanics on shotgun profiler tables
- experimental-design/multiple-testing - FDR control and multiplicity across taxa

---
name: bio-methylation-differential-cpg
description: Tests individual CpG sites for differential methylation (DMC/DMP) from bisulfite sequencing counts or array/continuous beta-value matrices. Covers the count-vs-continuous fork that dictates the model, beta-value vs M-value logit (Du 2010), beta-binomial overdispersion count models (DSS, methylKit, MOABS, RADMeth) for sequencing, limma moderated-t on M-values (eBayes trend/robust) for arrays, the bare-beta Welch t-test caveat, coverage-as-precision coupling, delta-beta effect size, BH-FDR with the neighboring-CpG dependence problem, EWAS genome-wide thresholds, and differential variability (DiffVar/iEVORA). Use when comparing per-CpG methylation between groups from WGBS/RRBS/targeted bisulfite or 450K/EPIC arrays, choosing a per-site test, or scanning for variance (not just mean) differences. For region-level aggregation see dmr-detection; for covariate/cell-fraction strategy and genomic inflation see ewas-design.
tool_type: mixed
primary_tool: limma
---

## Version Compatibility

Reference examples tested with: limma 3.58+, DSS 2.50+, methylKit 1.28+, missMethyl 1.36+, scipy 1.13+, statsmodels 0.14+, pandas 2.2+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

The DATA OBJECT is the version that matters most: sequencing yields integer (M, Cov) counts whose coverage is precision (a count model uses it); arrays yield a continuous beta with no coverage (a Gaussian model on M-values). The genome build of the calls (hg38 vs T2T-CHM13) and, for arrays, the platform (`450K` vs `EPIC`) fix the CpG universe and the genome-wide threshold. methylKit defaults (`overdispersion="none"`, `adjust="SLIM"`) silently change results; always confirm with `?calculateDiffMeth`.

# Per-CpG Differential Methylation Testing

**"Which single CpGs differ between my groups?"** -> First decide whether the data is sequencing COUNTS or a continuous array ratio, because that choice dictates the entire model - then test on M, report effect on beta, and gate on the intersection of FDR and |delta-beta|.
- R: `DSS::DMLtest()` on counts (sequencing); `limma::lmFit() |> eBayes(trend=TRUE, robust=TRUE)` on M-values (array/continuous)
- Python: `scipy.stats.ttest_ind(equal_var=False)` + `multipletests(method='fdr_bh')` - a continuous/array QUICK-LOOK only, never the headline sequencing test

Scope: the per-SITE test (DMC/DMP). Region-level aggregation (DSS callDMR, BSmooth, DMRcate) -> dmr-detection. Producing the (M, Cov) counts -> methylation-calling. Long-read MM/ML modBAM input -> long-read-sequencing/nanopore-methylation (pipe per-site counts back here). Covariate strategy, cell-fraction confounding, genomic inflation, replication design -> ewas-design.

## The Single Most Important Modern Insight -- The Right Test Is Dictated by the Data Object, and the Variance to Model Is Biological, Not Sampling Noise

The hardest lesson in the field is that the spread between biological replicates - not the coin-flip sampling at a single site - is the variance the test must capture. Three corollaries dictate the whole skill:

1. **Counts are not a beta.** Sequencing gives `(M, Cov)` per site per sample, and Cov IS precision: a beta of 0.80 from 80/100 reads is far more certain than 0.80 from 4/5. Collapsing to `beta = M/Cov` and running a t-test weights both sites equally and throws coverage away. For sequencing use a beta-binomial / overdispersion-corrected count model (DSS, methylKit `overdispersion="MN"`) that USES the coverage. Arrays genuinely have no counts - there a Gaussian model on M-values is correct.
2. **Binomial variance is not biological variance.** Fisher's exact (especially pooled across replicates) and uncorrected logistic regression assume the only randomness is binomial sampling at fixed depth. Two healthy individuals differ at a CpG far more than that, so these tests are anticonservative BY CONSTRUCTION - they hand back a long list of false positives that look exactly like findings. The whole job of DSS/MOABS/RADMeth and of methylKit's overdispersion option is to add the between-replicate (Beta) dispersion layer.
3. **Test on M, interpret on beta.** M-values (logit) are homoscedastic and well-calibrated; beta is bounded, heteroscedastic, and the only interpretable effect (delta-beta). The M-scale logFC is NOT a delta-beta and never maps linearly to one. The correct call is the INTERSECTION: `adj.P < cutoff AND |delta-beta| >= cutoff`.

Organize the analysis around defending these three, with the count-vs-continuous fork as the first decision - not around listing tests.

## Beta vs M-Value -- Why the Scale Matters

Beta = proportion methylated, range [0, 1], the unit of biological interpretation and of effect size (delta-beta). Its fatal property is heteroscedasticity: the variance of a proportion depends on its mean (maximal near 0.5, crushed toward 0 at the extremes where most genomic CpGs actually sit). A Gaussian linear model assumes constant variance, so a t-test/limma on raw beta is mis-calibrated, worst at the extremes. The M-value (Du 2010 *BMC Bioinformatics* 11:587), `log2(beta/(1-beta))`, is approximately homoscedastic and gives better-calibrated p-values, with the gap largest at high/low methylation. The division of labor: test on M, report delta-beta on beta. To avoid log(0) at beta in {0,1}, prefer computing M from intensities/counts as `log2((Meth+alpha)/(Unmeth+alpha))` (alpha ~ 1-100, never hits the boundary) over a symmetric offset on a precomputed beta.

## Tool Taxonomy

| Tool | Citation | Mechanism / role | When |
|------|----------|------------------|------|
| DSS | Feng 2014 *Nucleic Acids Res* 42:e69; Park & Wu 2016 *Bioinformatics* 32:1446 | beta-binomial, Bayesian dispersion shrinkage, Wald test | the count-based per-site default; few replicates; general designs |
| methylKit | Akalin 2012 *Genome Biol* 13:R87 | per-site logistic regression; `overdispersion="MN"` -> F-test | WGBS/RRBS; fast; SET overdispersion with replicates |
| MOABS (mcomp) | Sun 2014 *Genome Biol* 15:R38 | beta-binomial CDIF folding biological + statistical signal | CLI; depth-adjusted single metric |
| RADMeth | Dolzhenko & Smith 2014 *BMC Bioinformatics* 15:215 | beta-binomial regression, arbitrary multifactor design | CLI (methpipe); complex covariate models |
| limma | Ritchie 2015 *Nucleic Acids Res* 43:e47 | moderated-t on M-values, empirical-Bayes variance shrinkage | arrays (450K/EPIC) and any continuous matrix; small n |
| scipy Welch / Mann-Whitney | scipy docs | per-site continuous two-group test | array/continuous QUICK-LOOK only; NOT a count model |
| DiffVar (missMethyl) | Phipson & Oshlack 2014 *Genome Biol* 15:465 | Levene-style deviations + EB moderation (variance test) | scan for differentially VARIABLE CpGs alongside the mean |
| iEVORA | Teschendorff 2016 *Nat Commun* 7:10478 | Bartlett variance test + t re-ranking | field defects / rare stochastic outliers / risk prediction |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| WGBS/RRBS counts, replicates | DSS `DMLtest`, or methylKit `overdispersion="MN", test="F"` | beta-binomial uses coverage and models between-replicate dispersion |
| Sequencing, complex/multifactor design | DSS `DMLfit.multiFactor` or RADMeth | regression on counts with covariates |
| Unreplicated sequencing (n=1 vs n=1) | Fisher's exact on counts (exploratory) | no replicates means no biological variance to estimate; never pool replicates into this |
| 450K / EPIC array (or any continuous matrix) | limma moderated-t on M-values (`trend=TRUE, robust=TRUE`) | no counts exist; EB rescues small n; THE array workhorse |
| Quick continuous look, high uniform coverage | scipy Welch on M-values + BH | defensible shortcut at high depth; loses coverage-as-precision |
| Per-site, but biology is regional | per-site test, then -> dmr-detection | neighboring CpGs are correlated; aggregate for regional inference |
| Bulk tissue (blood etc.), any platform | add cell-fraction covariates to the design -> ewas-design | cell composition is the #1 EWAS confounder (Jaffe & Irizarry 2014) |
| Signal may be variance, not mean (cancer/aging) | DiffVar / iEVORA ON M-VALUES, alongside the mean scan | a variance change is invisible to a mean test |
| Long-read MM/ML per-site counts | -> long-read-sequencing/nanopore-methylation | calling/QC owned there; transfer (mod, valid) counts here |

When the design is contested (count model vs continuous, smooth vs not), verify current best practice against the installed tool's vignette rather than hard-coding one approach; count models are the rigorous default for sequencing, continuous tests a high-coverage shortcut, NOT a small-n shortcut.

## DSS Beta-Binomial on Counts (R, sequencing default)

**Goal:** Test each CpG for a mean methylation difference using a model that uses coverage as precision and shrinks the between-replicate dispersion.

**Approach:** Assemble per-sample (chr, pos, N=Cov, X=M) data frames into a BSseq object, run the Wald test per site with shrunken dispersion, then gate on FDR and an explicit effect-size floor.

```r
library(DSS)

# Each sample is a data.frame with columns chr, pos, N (total Cov), X (methylated M)
bs_obj <- makeBSseqData(list(c1, c2, c3, t1, t2, t3),
                        c('c1', 'c2', 'c3', 't1', 't2', 't3'))

# smoothing=FALSE keeps this a true per-SITE test; smoothing=TRUE borrows from
# neighbors (span 500 bp) and crosses into DMR territory -> dmr-detection
dml <- DMLtest(bs_obj, group1 = c('c1', 'c2', 'c3'), group2 = c('t1', 't2', 't3'), smoothing = FALSE)

# delta = effect-size floor on beta (default 0 applies NO gate); p.threshold is the FDR cut
dmc <- callDML(dml, delta = 0.1, p.threshold = 0.05)
# dml columns: mu1, mu2, diff (delta-beta on the beta scale), diff.se, stat, pval, fdr
```

## methylKit on Counts (R, the overdispersion trap)

**Goal:** Run a per-site logistic-regression test that accounts for between-replicate overdispersion (the default does not).

**Approach:** After uniting per-sample coverage objects, fit with overdispersion correction so the test becomes an F-test, then extract hyper/hypo sites with explicit effect and FDR floors.

```r
library(methylKit)

# meth is a united methylBase object (from methRead -> filterByCoverage -> unite)
# overdispersion="none" is the DEFAULT and over-calls under replication; set "MN" -> F-test
diff <- calculateDiffMeth(meth, overdispersion = 'MN', test = 'F', adjust = 'BH')
# adjust="SLIM" is the methylKit default, NOT BH; pass adjust="BH" to match other tools

# difference is in percentage points (25 = 25 points); qvalue is the FDR floor
dmc <- getMethylDiff(diff, difference = 25, qvalue = 0.01, type = 'all')
# meth.diff is a weighted-mean model difference, NOT mean(case_beta)-mean(ctrl_beta);
# recompute delta-beta from raw betas when comparing tools
```

## limma Moderated-t on M-Values (R, array/continuous default)

**Goal:** Identify DMPs from an array or continuous matrix at small n by borrowing variance across the ~10^5-10^6 probes.

**Approach:** Convert beta to M-values, fit per-probe linear models with EB moderation (trend+robust), extract BH-adjusted p-values, then attach delta-beta computed from the RAW betas.

```r
library(limma)

# M from intensities is cleaner; from a beta matrix use a boundary-safe transform
m_values <- log2((beta_matrix + 1e-3) / (1 - beta_matrix + 1e-3))

group <- factor(c(rep('case', 6), rep('ctrl', 6)))
design <- model.matrix(~ 0 + group)        # add cell-fraction/covariate columns here -> ewas-design
colnames(design) <- levels(group)
contrast_matrix <- makeContrasts(case - ctrl, levels = design)

fit <- lmFit(m_values, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)   # both default FALSE; set TRUE for methylation

res <- topTable(fit2, number = Inf, adjust.method = 'BH', sort.by = 'none')
# adjusted column is adj.P.Val (limma), NOT padj (DESeq2) or FDR; logFC is M-scale, NOT delta-beta
res$delta_beta <- rowMeans(beta_matrix[, group == 'case']) - rowMeans(beta_matrix[, group == 'ctrl'])
```

## Differential Variability -- Test the Variance, Not Just the Mean

A whole class of cancer/aging/field-defect signal lives in the SECOND moment: a CpG tight in controls (beta ~ 0.8) but scattered 0.3-0.95 in cases at the SAME mean is invisible to every mean test above. Run a differential-variability (DV) scan ALONGSIDE the mean scan; a CpG can be a DMP, a DVC, both, or neither.

**Goal:** Find CpGs whose spread (not mean) differs between groups, with FDR control robust to outliers.

**Approach:** On M-VALUES (the logit decouples variance from mean - see the boundary caveat below), fit Levene-style deviations with limma's EB moderation, then rank by adjusted p-value.

```r
library(missMethyl)

# Run on M-values: a variance difference on raw beta can be a pure mean-at-boundary artifact
fit <- varFit(m_values, design = design, coef = c(1, 2))   # ALWAYS pass coef (the group columns)
dvc <- topVar(fit, coef = 2, number = Inf)                 # coef must match; default is LAST column
# iEVORA (Bartlett variance + t re-ranking, Bartlett FDR < 0.001) is the field-defect alternative
```

Boundary caveat (the load-bearing DV trap): on beta in [0,1] variance is structurally tied to the mean (~ p(1-p)), so a mean shift from beta~0.95 toward ~0.6 mechanically RAISES variance and fakes a DV hit. Both DiffVar and iEVORA run on M-values for exactly this reason; always co-report the mean delta-beta next to any DV hit so a reader can judge whether the variance signal is independent of a boundary-driven mean move.

## Welch Quick-Look on Continuous Data (Python)

**Goal:** A fast continuous two-group per-site test for ARRAY/continuous matrices (or high uniform-coverage sequencing where coverage loss is accepted) - explicitly not the sequencing headline.

**Approach:** Test on M-values per CpG with Welch (unequal variance), then apply BH FDR; report delta-beta from raw betas.

```python
import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

# m_case / m_ctrl: M-value matrices (rows CpGs, cols samples). For SEQUENCING counts prefer DSS.
_, pvalues = ttest_ind(m_case, m_ctrl, axis=1, equal_var=False, nan_policy='omit')  # Welch; scipy default is Student's
reject, padj, _, _ = multipletests(pvalues, method='fdr_bh')  # default is 'hs' (Holm-Sidak), must set fdr_bh
delta_beta = beta_case.mean(axis=1) - beta_ctrl.mean(axis=1)  # effect on the beta scale
```

## Multiple Testing and the Dependence Problem

BH-FDR is the default at both platforms; Bonferroni leaves almost nothing at 850k EPIC probes or 28M+ WGBS CpGs. Two caveats: (1) BH assumes independence (or positive dependence), but neighboring CpGs are strongly spatially correlated, so per-site BH on methylation is conservative-but-not-exact, and a lone significant CpG flanked by null neighbors is suspect - this regional dependence is precisely why region-level methods exist (-> dmr-detection); do not try to fix it inside per-site BH. (2) Large consortium array EWAS often use a fixed genome-wide threshold for comparability instead of BH: ~2.4e-7 experiment-wide for 450K (Saffari 2018), and P < 9e-8 (~8.6e-9 genome-wide) for EPIC (Mansell 2019). WGBS has no single accepted constant - BH or region-level FDR dominate.

## Per-Method Failure Modes

### Bare-beta t-test on sequencing counts
**Trigger:** computing `beta = M/Cov` from bisulfite counts and running a t-test. **Mechanism:** discards coverage (an 8-read and an 800-read site weigh equally) and, on raw beta, fights heteroscedasticity. **Symptom:** noisy low-coverage sites masquerade as confident hits; poor replication. **Fix:** for sequencing use DSS or methylKit `overdispersion="MN"` on counts; if forced continuous, at least test on M-values at high depth.

### methylKit left at overdispersion="none"
**Trigger:** `calculateDiffMeth(meth)` with replicates and no overdispersion argument. **Mechanism:** the default does NO overdispersion correction - a plain logistic LRT that assumes binomial-only variance. **Symptom:** inflated significant-CpG count; anticonservative p-values. **Fix:** `overdispersion="MN", test="F"`; pass `adjust="BH"` (default is SLIM).

### Pooling replicates for Fisher's exact
**Trigger:** summing M and U across replicates into one super-sample per group to "have enough counts". **Mechanism:** collapses biological variance - treats N mice as one giant mouse. **Symptom:** wildly anticonservative p-values. **Fix:** a replicate-aware count model; Fisher only for n=1 vs n=1, reported as exploratory.

### Reporting the M-scale logFC as delta-beta
**Trigger:** quoting limma `logFC` or methylKit `meth.diff` as the methylation-percentage change. **Mechanism:** the logit is steep at 0.5 and flat at the ends, so the same logFC is a large beta-change mid-range and a tiny one at the extremes. **Symptom:** overstated effect sizes near 0/1. **Fix:** recompute delta-beta from raw betas for reporting, always.

### Significance without effect size (and the reverse)
**Trigger:** ranking by FDR alone at large n, or by delta-beta alone at small n. **Mechanism:** at 28M CpGs a 2-point delta-beta within noise clears FDR; a big delta from 3 noisy samples is not a finding. **Symptom:** non-reproducible top hits. **Fix:** gate on the intersection `adj.P < cutoff AND |delta-beta| >= cutoff`.

### Winner's curse on discovery effect sizes
**Trigger:** quoting the top hits' discovery delta-betas as the true effect. **Mechanism:** thresholding selects sites where noise pushed the estimate up, so reported magnitudes are upward-biased. **Symptom:** replication cohorts show attenuated effects; replications powered on the inflated delta-beta are underpowered. **Fix:** flag discovery |delta-beta| as an upper bound; estimate the honest effect from independent replication (Palmer & Pe'er 2017).

### Ignoring cell composition in bulk tissue
**Trigger:** an EWAS on whole blood / bulk tissue with no cell-fraction covariates. **Mechanism:** the top "DMPs" are often shifts in cell-type proportion, not within-cell methylation (Jaffe & Irizarry 2014). **Symptom:** hits that fail to replicate across cohorts. **Fix:** estimate cell fractions and add them to the design matrix -> ewas-design.

### Differential-variability hit that is a mean artifact
**Trigger:** a DV test on raw beta values. **Mechanism:** beta variance is tied to the mean, so a mean move toward 0.5 fakes higher variance. **Symptom:** DV hits that co-occur with large mean shifts toward 0.5. **Fix:** run DiffVar/iEVORA on M-values and co-report the mean delta-beta.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Min coverage 10x in EVERY sample | field standard | below it a single-CpG beta is granular and noisy; floor must hold in all compared samples |
| Upper cap at 99.9th percentile Cov | field standard | drops PCR-duplicate pileups / collapsed-repeat mapping artifacts |
| WGBS 5-10x / RRBS 10x / targeted 30-100x | assay convention | smoothing tolerates 5x; targeted expects deep, even depth |
| delta-beta floor 0.10 / 0.20 / 0.30 | convention | 0.10 EWAS discovery (diluted by cell mixture), 0.20 general, 0.30 cancer-vs-normal |
| methylKit getMethylDiff difference=25, qvalue=0.01 | Akalin 2012 *Genome Biol* 13:R87 | tool defaults; 25 percentage points, FDR 0.01 |
| DSS callDML delta default 0 (set it) | DSS docs | delta=0 applies NO effect-size gate; set delta=0.1 |
| BH-FDR, not Bonferroni | field standard | Bonferroni leaves almost nothing at 10^5-10^7 tests |
| 450K ~2.4e-7; EPIC P<9e-8 | Saffari 2018; Mansell 2019 | fixed array EWAS thresholds for cross-study comparability |
| Fisher's exact: n=1 vs n=1 only | mechanism | no replicates means no biological variance; never pool replicates |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Long list of low-coverage false positives | bare-beta t-test on sequencing counts | DSS / methylKit `overdispersion="MN"` on counts |
| methylKit returns too many DMCs | left at `overdispersion="none"` | set `overdispersion="MN", test="F"` |
| methylKit q-values disagree with other tools | default `adjust="SLIM"`, not BH | pass `adjust="BH"` |
| Effect sizes overstated near 0/1 | reported M-scale logFC as delta-beta | recompute delta-beta from raw betas |
| Almost nothing significant genome-wide | Bonferroni at millions of tests | BH-FDR (WGBS) or EWAS thresholds (array) |
| `padj`/`$FDR` column not found in limma | wrong column name | limma column is `adj.P.Val` |
| All p-values NaN in Python | `multipletests` default `method='hs'` or wrong axis | pass `method='fdr_bh'`; test on M-values, `axis=1` |
| EWAS hits do not replicate | cell composition / winner's curse | cell-fraction covariates; treat discovery effects as upper bounds |

## References

- Du P, Zhang X, Huang C-C, Jafari N, Kibbe WA, Hou L, Lin SM. 2010. Comparison of Beta-value and M-value methods for quantifying methylation levels by microarray analysis. *BMC Bioinformatics* 11:587.
- Feng H, Conneely KN, Wu H. 2014. A Bayesian hierarchical model to detect differentially methylated loci from single nucleotide resolution sequencing data. *Nucleic Acids Res* 42:e69.
- Park Y, Wu H. 2016. Differential methylation analysis for BS-seq data under general experimental design. *Bioinformatics* 32:1446-1453.
- Akalin A, Kormaksson M, Li S, Garrett-Bakelman FE, Figueroa ME, Melnick A, Mason CE. 2012. methylKit: a comprehensive R package for the analysis of genome-wide DNA methylation profiles. *Genome Biol* 13:R87.
- Sun D, Xi Y, Rodriguez B, Park HJ, Tong P, Meong M, Goodell MA, Li W. 2014. MOABS: model based analysis of bisulfite sequencing data. *Genome Biol* 15:R38.
- Dolzhenko E, Smith AD. 2014. Using beta-binomial regression for high-precision differential methylation analysis in multifactor whole-genome bisulfite sequencing experiments. *BMC Bioinformatics* 15:215.
- Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK. 2015. limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Res* 43:e47.
- Phipson B, Oshlack A. 2014. DiffVar: a new method for detecting differential variability with application to methylation in cancer and aging. *Genome Biol* 15:465.
- Teschendorff AE, Gao Y, Jones A, Ruebner M, Beckmann MW, Wachter DL, Fasching PA, Widschwendter M. 2016. DNA methylation outliers in normal breast tissue identify field defects that are enriched in cancer. *Nat Commun* 7:10478.
- Jaffe AE, Irizarry RA. 2014. Accounting for cellular heterogeneity is critical in epigenome-wide association studies. *Genome Biol* 15:R31.
- Saffari A, Silver MJ, Zavattari P, Moi L, Columbano A, Meaburn EL, Dudbridge F. 2018. Estimation of a significance threshold for epigenome-wide association studies. *Genet Epidemiol* 42:20-33.
- Mansell G, Gorrie-Stone TJ, Bao Y, Kumari M, Schalkwyk LS, Mill J, Hannon E. 2019. Guidance for DNA methylation studies: statistical insights from the Illumina EPIC array. *BMC Genomics* 20:366.
- Palmer C, Pe'er I. 2017. Statistical correction of the Winner's Curse explains replication variability in quantitative trait genome-wide association studies. *PLoS Genet* 13:e1006916.

## Related Skills

- methylation-calling - Produces the (M, coverage) counts tested here
- methylkit-analysis - methylKit object model and calculateDiffMeth mechanics
- dmr-detection - Region-level aggregation downstream of per-site testing
- cell-type-deconvolution - Cell-fraction covariates (the dominant bulk-tissue confounder)
- ewas-design - Covariate strategy, genomic inflation, and genome-wide thresholds
- experimental-design/multiple-testing - FDR/FWER theory behind the corrections applied here
- long-read-sequencing/nanopore-methylation - Long-read MM/ML calling; pipe per-site counts here for count-based statistics
- differential-expression/deseq2-basics - Analogous dispersion-shrinkage / empirical-Bayes machinery
- workflows/methylation-pipeline - End-to-end bisulfite pipeline

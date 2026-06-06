---
name: bio-metabolomics-statistical-analysis
description: Decision-grade statistical analysis for metabolomics intensity tables. Covers transformation and scaling (Pareto vs unit-variance as a hidden hypothesis), unsupervised structure (PCA/HCA for QC), permutation-validated PLS-DA/OPLS-DA (R2 vs Q2, double CV, VIP as heuristic), univariate testing (Welch/Mann-Whitney/ANOVA/LMM with covariate adjustment), and dependence-aware multiple testing. Use when testing which metabolites differ, building or validating a discriminant model, choosing a scaling, or correcting many correlated tests. For sample-wise normalization/drift correction see metabolomics/normalization-qc; for ML classifiers and selection-inside-CV leakage see machine-learning/biomarker-discovery and machine-learning/model-validation; for pathway interpretation see metabolomics/pathway-mapping; for design/power/multiplicity regime see experimental-design/multiple-testing.
tool_type: mixed
primary_tool: ropls
---

## Version Compatibility

Reference examples tested with: ropls 1.34+, scipy 1.12+, statsmodels 0.14+, numpy 1.26+, pandas 2.1+, matplotlib 3.8+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Metabolomics Statistical Analysis

**"Tell me which metabolites separate my groups"** -> run an honest univariate test with dependence-aware FDR AND a permutation-validated multivariate model, then reconcile the two.
- R: `ropls::opls()` (PCA/PLS-DA/OPLS-DA + permutation), `wilcox.test()`/`lm()`/`lme4::lmer()`, `p.adjust(method='BH')`
- Python: `scipy.stats.ttest_ind(equal_var=False)`/`mannwhitneyu`, `statsmodels` `multipletests(method='fdr_bh')`, `sklearn.cross_decomposition.PLSRegression` + `permutation_test_score`

## The Single Most Important Modern Insight -- A Score Plot Is a Hypothesis, Not a Result

In metabolomics the regime is p >> n (hundreds-to-thousands of features, tens of samples) with strongly correlated features. In that regime any binary labelling of n points in >= n-1 dimensions is linearly separable with probability 1, so PLS-DA and OPLS-DA produce a clean two-cluster score plot even for randomly assigned labels. A beautiful score plot is the generic output of the algorithm and carries essentially zero information. Only a cross-validated Q2 benchmarked against a permutation null distinguishes signal from geometry (Westerhuis 2008; Ruiz-Perez 2020). Three corollaries reorganize the whole skill: (1) R2 is no evidence (it can be driven to 1 by adding components); only permutation-validated Q2 licenses a claim. (2) Scaling is a hidden hypothesis -- variance-driven methods weight a feature by the variance it is allowed to contribute, so Pareto vs unit-variance hands back a different VIP list and a different biological story (van den Berg 2006). (3) Features are not independent (pathways, adducts, isotopologues), so naive BH-independence is violated and one real signal lights up its whole correlated cluster.

## Transformation and Scaling -- the Decision That Changes Conclusions

Transformation (nonlinear, per value: corrects heteroscedastic multiplicative MS noise and skew) and scaling (linear, per feature: sets relative weight) are distinct. Mean-centering is the universal first step. Choosing not to scale is the strongest prior of all -- it lets the most abundant metabolite drive PC1.

| Method | Per feature j | Effect | Use when |
|--------|---------------|--------|----------|
| Centering | subtract mean | offsets removed, variance unchanged | always (prerequisite for all below) |
| Auto / unit-variance / "standard" | center, / SD_j | every metabolite equal weight | a priori all metabolites equally important; classic default -- but inflates near-LOD noise |
| Pareto | center, / sqrt(SD_j) | between raw and UV | de facto metabolomics/NMR/OPLS-DA default; curbs dominance with less noise inflation than UV |
| Range | center, / (max-min) | abundance dependence removed | clean data, few outliers (outlier-sensitive) |
| Vast | UV x (mean/SD) | up-weights low-CV stable features | focus on robust/reproducible features with prior class info |
| Level | center, / mean_j | relative (% change) | response as relative change (mean noisy at low abundance) |
| Log / log10 | log(x) | multiplicative -> additive | concentrations spanning orders of magnitude (undefined at 0) |
| glog | linear near 0, log for large x | variance-stabilizing | data with zeros / near-LOD values; preferred over plain log (needs transition param) |
| Power (sqrt, cube-root) | x^(1/2) | mild stabilization | mild skew with zeros present |

van den Berg 2006: on real data autoscale and range recovered biologically meaningful loadings; Pareto is the pragmatic middle. Decision rule: transform first (if heteroscedastic -- usually yes for MS), then center, then scale; run at least Pareto AND UV, and if the top-VIP list or conclusion flips, the result is scaling-fragile and must be tempered.

## Decision Tree by Scenario

| Goal / situation | Do | Why |
|------------------|----|----|
| First look, QC, batch/outlier check | PCA on scaled data; color scores by batch/injection order; Hotelling T2 ellipse | Unsupervised -> cannot overfit the grouping; pooled-QC samples must cluster tightly in the center, else analytical variance dominates |
| Which single metabolites differ (2 groups) | Welch t-test (post-transform) or Mann-Whitney; BH FDR; report fold change + CI | Interpretable per-feature effect; FDR-controlled; effect size mandatory in p>>n |
| 2 groups, paired/pre-post | Paired t-test or Wilcoxon signed-rank | Discards within-subject pairing if analyzed unpaired -> underpowered |
| >2 groups | One-way ANOVA (+Tukey) or Kruskal-Wallis (+Dunn) | Match normality assumption |
| Longitudinal / repeated measures | Linear mixed model (random intercept/slope per subject) | Handles unbalanced timepoints, missingness, within-subject correlation |
| Covariate adjustment (age/sex/BMI/batch) | Per-feature linear model `y ~ group + covars` | Human metabolome is dominated by age/sex/BMI -- unadjusted, they masquerade as case/control signal (Thevenot 2015) |
| A discriminant / predictive model | PLS-DA (`orthoI=0`) or OPLS-DA (`predI=1, orthoI=NA`) + permutation + double CV | Supervised; demands full validation (see checklist) |
| Built-in feature selection | sparse PLS-DA `splsda()` (mixOmics) with `tune.splsda` | Selection must be inside CV -> hand off to machine-learning/biomarker-discovery |
| Rank discriminant features | VIP from a permutation-validated model only; corroborate with univariate FDR | VIP > 1 is a heuristic, not a test (see failure modes) |
| Confirm a biomarker | Independent validation cohort | Internal CV does not correct overfitting/forking-paths; discovery performance overestimates external |

## Scaling + PCA

**Goal:** Get the honest unsupervised first look that cannot chase the labels, with QC as the primary data-quality readout.

**Approach:** Transform if heteroscedastic, then PCA with an explicit scaling; inspect QC clustering, batch coloring, and Hotelling T2.

```r
library(ropls)
# scaleC default is "standard" (unit-variance/autoscale), NOT Pareto -- set explicitly
pca <- opls(t(feature_matrix), scaleC = 'pareto', fig.pdfC = 'none', info.txtC = 'none')
scores <- getScoreMN(pca)               # samples x components
getSummaryDF(pca)                       # R2X(cum) per component
# Tight pooled-QC clustering in the center = trustworthy run; QC scatter = analytical variance dominates
```

## Permutation-Validated PLS-DA / OPLS-DA

**Goal:** Decide whether group separation is real, not a geometry artifact, before reading any VIP or S-plot.

**Approach:** Fit with an explicit scaling, raise `permI` far above the default of 20, and read `pQ2`/`pR2Y` -- a model whose true Q2 sits inside the permutation cloud is indistinguishable from chance.

```r
library(ropls)
group <- factor(sample_info$group)
# OPLS-DA: 1 predictive + auto orthogonal; permI default 20 is too few for a reliable pQ2 -> >=1000
oplsda <- opls(t(feature_matrix), group, predI = 1, orthoI = NA,
               scaleC = 'pareto', permI = 1000, crossvalI = 7,
               fig.pdfC = 'none', info.txtC = 'none')
summ <- getSummaryDF(oplsda)            # R2X(cum), R2Y(cum), Q2(cum), pre, ort, pR2Y, pQ2
vip_pred <- getVipVn(oplsda)            # predictive VIP (Galindo-Prieto 2014); orthoL=TRUE for orthogonal
# Claim is licensed only if Q2 high AND pQ2 small. R2Y alone proves nothing.
```

PLS-DA is `orthoI = 0`. OPLS-DA has identical predictive power to PLS-DA -- it is a coordinate rotation, not a better model; the orthogonal block often encodes a confounder (inspect what correlates with it). DQ2 (Westerhuis 2008b) is the discriminant-appropriate figure of merit when Q2 penalizes correct-side over-predictions.

## PLS-DA / OPLS-DA Validation Checklist

1. Report the transformation + scaling used (it changes the loadings, VIPs, and story).
2. Report R2X, R2Y, Q2 and the number of predictive + orthogonal components.
3. Choose the number of components inside CV, not by eye on the training fit.
4. Permutation test (>= 1000) of the full pipeline -> permutation p for Q2 (and R2Y). Permute every step that touched the labels.
5. For honest generalization error use double (cross-model) CV or an untouched test set; single CV that also tuned the model is optimistic.
6. Independent validation cohort for any biomarker claim.
7. Read VIP / S-plot only from a validated model; corroborate with univariate FDR + effect size; report ranking stability across resamples.
8. Put a PCA score plot beside the PLS-DA one -- separation only under supervision is the artifact signature.

## Univariate Testing + Correct FDR

**Goal:** Produce an interpretable, FDR-controlled per-metabolite answer with effect sizes.

**Approach:** Match the test to the design, compute log2 fold change as a difference of group means on transformed data, then apply BH explicitly (defaults are not BH in either language).

```python
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

logged = np.log2(intensities.replace(0, np.nan))   # transform before testing
pvals, lfc = [], []
for feat in logged.index:
    a = logged.loc[feat, case].dropna().values
    b = logged.loc[feat, ctrl].dropna().values
    if len(a) >= 3 and len(b) >= 3:
        pvals.append(ttest_ind(a, b, equal_var=False)[1])   # Welch: scipy defaults to Student
        lfc.append(a.mean() - b.mean())                     # geometric-mean ratio on log scale
    else:
        pvals.append(np.nan); lfc.append(np.nan)
res = pd.DataFrame({'feature': logged.index, 'log2fc': lfc, 'pval': pvals}).dropna(subset=['pval'])
# statsmodels default is 'hs' (Holm-Sidak); R p.adjust default is 'holm' -- ALWAYS pass BH explicitly
res['padj'] = multipletests(res['pval'], method='fdr_bh')[1]
```

BH controls FDR under independence and PRDS; positively-correlated metabolomics features roughly satisfy PRDS, so BH is valid but conservative -- but closure-induced negative correlations (after total-area/PQN normalization) fall outside the clean case, where a permutation FDR sidesteps the dependence assumptions. The effective number of independent tests is far below the feature count (one compound = many adducts/isotopologues/fragments); use an effective-number-of-tests correction (Peluso 2021) rather than Bonferroni-on-features, and collapse features to compounds before counting "how many metabolites changed."

## Volcano Plot

**Goal:** Show significance and magnitude together for all features.

**Approach:** Plot log2 fold change vs -log10(p), with the FDR cutoff annotated (raw p on the axis is fine only if the FDR line is drawn).

```python
import matplotlib.pyplot as plt
hit = (res['padj'] < 0.05) & (res['log2fc'].abs() > 1)   # 2-fold + FDR 5%
plt.scatter(res['log2fc'], -np.log10(res['pval']), c=np.where(hit, 'firebrick', 'gray'), s=12, alpha=0.6)
plt.axhline(-np.log10(0.05), ls='--'); plt.axvline(1, ls='--'); plt.axvline(-1, ls='--')
plt.xlabel('log2 fold change'); plt.ylabel('-log10(p)')
```

## Per-Method Failure Modes

### Noise separation (the cardinal sin)
- **Trigger:** Reporting a PLS-DA/OPLS-DA score plot as evidence of a group difference.
- **Mechanism:** In p>>n any labelling is linearly separable; the algorithm always finds a covariance-maximizing direction, even for random labels.
- **Symptom:** Clean two-cluster score plot, high R2Y, but Q2 low/negative or inside the permutation cloud; PCA shows no separation.
- **Fix:** Permutation test (>=1000) of the full pipeline; require small pQ2; put the PCA plot beside it.

### VIP misuse
- **Trigger:** Selecting biomarkers by VIP > 1 from a single model fit.
- **Mechanism:** VIPs are normalized so the mean squared VIP = 1 -- roughly half the features exceed 1 by construction; there is no null, no p-value, and the ranking is unstable under resampling in p>>n.
- **Symptom:** Top-20 VIP list reshuffles when the model is re-bootstrapped; VIP-only hits fail to replicate.
- **Fix:** Use VIP only from a permutation-validated model; require univariate FDR + effect-size concordance and resampling stability; use the OPLS-specific VIP so a high orthogonal-block VIP (the confounder) is not credited to disease.

### Naive FDR under correlation
- **Trigger:** BH or Bonferroni applied as if the features were independent.
- **Mechanism:** Pathway co-regulation plus adducts/isotopologues/fragments make features strongly correlated; one signal lights up its whole cluster, and closure (after sample-wise normalization) injects negative correlations.
- **Symptom:** A "200 significant metabolites" list that encodes a handful of independent signals; over-conservative threshold from Bonferroni-on-features.
- **Fix:** Effective-number-of-tests or permutation FDR (Peluso 2021); collapse features to compounds before counting hits; report independent-signal counts.

### Log with zeros / detection-rate confound
- **Trigger:** Half-min (or zero) imputation followed by log, especially when detection rate differs between groups.
- **Mechanism:** "Missing" is left-censored (MNAR); a constant imputed at the LOD then logged spikes the censored region, and a detection-rate difference masquerades as a concentration difference.
- **Symptom:** Fake bimodality; a low-abundance "hit" that is really a difference in how often the metabolite was detected.
- **Fix:** Report per-group detection rates with any low-abundance hit; prefer glog or a left-censored imputer (QRILC/GSimp) over impute-constant-then-log when detection differs (see metabolomics/normalization-qc).

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Q2 > 0.5 "good" | Triba 2015 (heuristic) | Predictive ability rule-of-thumb; not a hard cutoff -- many published models report Q2 < 0.5; report the value, not a verdict |
| permI >= 1000 | Szymanska 2012 | Q2/DQ2 null distributions are skewed; the ropls default of 20 estimates only the granularity of the grid, not a usable pQ2 |
| pQ2 < 0.05 | Westerhuis 2008 | Fraction of permuted models with Q2 >= true Q2; the actual evidence the separation is real |
| crossvalI = 7 | ropls default | 7-fold CV; for very small n LOO is common but optimistic |
| VIP > 1 | Galindo-Prieto 2014 | Above-average contributor; a ranking heuristic with no error control -- never a standalone selector |
| BH FDR < 0.05 | Benjamini-Hochberg | Expected false-positive proportion among rejections; the metabolomics discovery default |
| \|log2FC\| > 1 | convention | 2-fold; effect-size gate orthogonal to the p-value, mandatory in p>>n |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Model "significant" yet noise | `permI = 20` (ropls default) | Set `permI >= 1000`; read `pQ2`/`pR2Y` from `getSummaryDF` |
| Wrong scaling shipped silently | `scaleC` default is `"standard"` (UV), not Pareto | Set `scaleC = 'pareto'` (or the intended scaling) explicitly; report it |
| PLS-DA vs OPLS-DA "function not found" | type is set by `orthoI`, not a separate function | `orthoI = 0` -> PLS; `orthoI = NA` -> OPLS; `predI = 1` for 2-class |
| FDR is actually Holm | R `p.adjust` default is `'holm'` (FWER) | Pass `method = 'BH'` |
| FDR is actually Holm-Sidak | statsmodels `multipletests` default is `'hs'` | Pass `method = 'fdr_bh'` |
| Student instead of Welch | scipy `ttest_ind` default `equal_var=True` | Set `equal_var=False` (group variances differ, esp. near LOD) |
| Reversed/unstable fold change | `log2(mean_ratio)` uses arithmetic means | Difference of log-means (geometric-mean ratio), consistent with limma/DESeq2 |
| Optimistic CV error | feature selection done before CV | Re-fit selection inside every fold; see machine-learning/model-validation |
| getVipVn gives orthogonal importance | `orthoL = TRUE` returns orthogonal VIP | Use default (predictive VIP) for discriminant ranking |

## References

- van den Berg RA, Hoefsloot HCJ, Westerhuis JA, Smilde AK, van der Werf MJ. 2006. Centering, scaling, and transformations: improving the biological information content of metabolomics data. *BMC Genomics* 7:142.
- Westerhuis JA, Hoefsloot HCJ, Smit S, Vis DJ, Smilde AK, et al. 2008. Assessment of PLSDA cross validation. *Metabolomics* 4:81-89.
- Westerhuis JA, van Velzen EJJ, Hoefsloot HCJ, Smilde AK. 2008. Discriminant Q2 (DQ2) for improved discrimination in PLSDA models. *Metabolomics* 4:293-296.
- Saccenti E, Hoefsloot HCJ, Smilde AK, Westerhuis JA, Hendriks MMWB. 2014. Reflections on univariate and multivariate analysis of metabolomics data. *Metabolomics* 10:361-374.
- Broadhurst DI, Kell DB. 2006. Statistical strategies for avoiding false discoveries in metabolomics and related experiments. *Metabolomics* 2:171-196.
- Thevenot EA, Roux A, Xu Y, Ezan E, Junot C. 2015. Analysis of the human adult urinary metabolome variations with age, body mass index, and gender by implementing a comprehensive workflow for univariate and OPLS statistical analyses. *J Proteome Res* 14:3322-3335.
- Triba MN, Le Moyec L, Amathieu R, Goossens C, Bouchemal N, et al. 2015. PLS/OPLS models in metabolomics: the impact of permutation of dataset rows on the K-fold cross-validation quality parameters. *Mol BioSyst* 11:13-19.
- Szymanska E, Saccenti E, Smilde AK, Westerhuis JA. 2012. Double-check: validation of diagnostic statistics for PLS-DA models in metabolomics studies. *Metabolomics* 8(Suppl 1):3-16.
- Galindo-Prieto B, Eriksson L, Trygg J. 2014. Variable influence on projection (VIP) for orthogonal projections to latent structures (OPLS). *J Chemometr* 28:623-632.
- Ruiz-Perez D, Guan H, Madhivanan P, Mathee K, Narasimhan G. 2020. So you think you can PLS-DA? *BMC Bioinformatics* 21(Suppl 1):2.
- Peluso A, Glen R, Ebbels TMD. 2021. Multiple-testing correction in metabolome-wide association studies. *BMC Bioinformatics* 22:67.
- Storey JD, Tibshirani R. 2003. Statistical significance for genomewide studies. *Proc Natl Acad Sci USA* 100:9440-9445.

## Related Skills

- metabolomics/normalization-qc - Sample-wise normalization, drift/batch correction, missing-value imputation upstream of testing
- metabolomics/pathway-mapping - Functional interpretation of differential metabolites
- machine-learning/biomarker-discovery - Feature selection inside CV, stability, minimal-optimal vs all-relevant
- machine-learning/model-validation - Leakage taxonomy, nested CV, calibration vs discrimination
- experimental-design/multiple-testing - FDR vs FWER regime, discovery vs confirmatory
- data-visualization/volcano-and-ma-plots - Volcano plot recipes

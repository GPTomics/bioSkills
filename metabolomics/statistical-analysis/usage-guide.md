# Statistical Analysis Usage Guide

## Overview

Statistical analysis decides which metabolites differ between conditions and whether a discriminant model is real. The two endemic failure modes it guards against are overfitting (in p >> n, supervised models separate pure noise -- a score plot proves nothing without permutation-validated Q2) and false-discovery inflation (thousands of correlated tests on small n yield long "significant" lists that do not replicate). It covers transformation/scaling as an explicit modeling choice, unsupervised structure (PCA/HCA for QC), validated PLS-DA/OPLS-DA, univariate testing with covariate adjustment, and dependence-aware multiple testing.

## Prerequisites

```bash
# Python
pip install scipy statsmodels numpy pandas matplotlib scikit-learn
```

```r
# R
BiocManager::install("ropls")
install.packages(c("mixOmics", "lme4", "qvalue"))
```

Conceptual prerequisites: decide the scaling (Pareto vs unit-variance is a hypothesis, not a default), know whether the normalization imposed closure (total-area/PQN), distinguish technical from biological zeros, and identify the experimental unit (the subject, not the injection) before any group test.

## Quick Start

Tell your AI agent what you want to do:
- "Run PCA with Pareto scaling and check whether my pooled QC samples cluster tightly"
- "Fit an OPLS-DA model and permutation-test it -- I only trust it if pQ2 is small"
- "Run Welch t-tests between case and control with explicit BH FDR and fold changes"
- "Fit a linear mixed model per metabolite for my longitudinal design"
- "Show me whether my top hits survive switching from Pareto to unit-variance scaling"

## Example Prompts

### Scaling and Unsupervised Structure
> "Transform and Pareto-scale my feature table, run PCA, and color the scores by batch and injection order to check for drift."
> "Re-run the PCA with unit-variance scaling and tell me whether the top loadings change."

### Validated Multivariate Models
> "Build an OPLS-DA model between disease and control, run 1000 permutations, and report R2X, R2Y, Q2, and pQ2."
> "Put a PCA score plot next to my PLS-DA plot so I can see whether the separation is real or supervised-only."
> "Extract VIP scores but only after the model passes the permutation test, and cross-check them against my univariate hits."

### Univariate Testing and FDR
> "Run Welch t-tests in Python with Benjamini-Hochberg correction and report log2 fold changes."
> "Fit a per-metabolite linear model adjusting for age, sex, and BMI."
> "My features are heavily correlated -- use an effective-number-of-tests correction instead of Bonferroni."

### Reproducibility
> "Check whether my discriminant result survives changing the scaling and the imputation method."
> "I have a validation cohort -- estimate external performance, not just cross-validated discovery performance."

## What the Agent Will Do

1. Confirm the upstream state: normalization/closure, zero handling, and the experimental unit.
2. Choose and report a transformation + scaling, and offer to run a second scaling as a robustness check.
3. Run PCA for QC and structure (pooled-QC clustering, batch/outlier inspection) before any supervised model.
4. Match the univariate test to the design (Welch/Mann-Whitney/ANOVA/paired/LMM), with covariate adjustment where needed.
5. Apply BH FDR explicitly (defaults are not BH in R or statsmodels) and report fold change + CI; use a dependence-aware correction when features are strongly correlated.
6. For any discriminant model, fit PLS-DA/OPLS-DA, raise permutations to >= 1000, and report pQ2 plus the validation checklist.
7. Read VIP/S-plot only from a validated model and reconcile multivariate VIPs with univariate FDR.
8. Generate a volcano plot and a results table; flag scaling-fragile or detection-rate-confounded hits.

## Tips
- There is no scaling-free analysis: report the scaling, and if conclusions flip between Pareto and UV, the result is a property of the prior, not the biology.
- `ropls::opls()` defaults `scaleC = "standard"` (unit-variance), not Pareto -- set it explicitly.
- Raise `permI` well above the ropls default of 20 (>= 1000) or the permutation p-value just reports the grid granularity.
- R `p.adjust` defaults to Holm (FWER) and statsmodels `multipletests` to Holm-Sidak -- always pass BH explicitly.
- In Python set `equal_var=False` in `scipy.stats.ttest_ind()` for Welch; group variances differ, especially near the LOD.
- Compute fold change as a difference of means on the log scale (geometric-mean ratio), not `log2(mean_ratio)`.
- VIP > 1 is a heuristic with no error control; corroborate with univariate FDR + effect size and resampling stability.
- Internal cross-validation does not substitute for an independent cohort -- discovery performance overestimates external.
- Report per-group detection rates with any low-abundance hit; a detection-rate difference can masquerade as a concentration difference.

## Related Skills

- metabolomics/normalization-qc - Sample-wise normalization, drift/batch correction, missing-value imputation upstream of testing
- metabolomics/pathway-mapping - Functional interpretation of differential metabolites
- machine-learning/biomarker-discovery - Feature selection inside CV, stability, minimal-optimal vs all-relevant
- machine-learning/model-validation - Leakage taxonomy, nested CV, calibration vs discrimination
- experimental-design/multiple-testing - FDR vs FWER regime, discovery vs confirmatory
- data-visualization/volcano-and-ma-plots - Volcano plot recipes

## References

- van den Berg RA, et al. 2006. Centering, scaling, and transformations. *BMC Genomics* 7:142.
- Westerhuis JA, et al. 2008. Assessment of PLSDA cross validation. *Metabolomics* 4:81-89.
- Thevenot EA, et al. 2015. ropls / urinary metabolome workflow. *J Proteome Res* 14:3322-3335.
- Szymanska E, et al. 2012. Double-check: validation of diagnostic statistics for PLS-DA. *Metabolomics* 8(Suppl 1):3-16.
- Peluso A, Glen R, Ebbels TMD. 2021. Multiple-testing correction in metabolome-wide association studies. *BMC Bioinformatics* 22:67.

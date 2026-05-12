# Mediation Analysis - Usage Guide

## Overview

Decomposes the total effect of an exposure (genotype, treatment, environmental factor) on an outcome into direct and indirect paths through one or more mediators. Covers single-mediator observational mediation (`mediation::mediate`), 4-way decomposition with exposure-mediator interaction (CMAverse), high-dimensional EWAS / transcriptome-wide mediator screening (HIMA / HIMA2 / BAMA), MR-based mediation when sequential ignorability is implausible (two-step MR, MVMR-mediation), and doubly-robust double-ML mediation (`causalweight::medDML`). Sequential ignorability is fundamentally untestable, so every reported result is paired with a sensitivity analysis (Imai rho or mediational E-value).

## Prerequisites

```r
# CRAN
install.packages(c('mediation','HIMA','bama','causalweight','EValue'))

# GitHub (NOT on CRAN)
remotes::install_github('BS1125/CMAverse')
remotes::install_github('WSpiller/MVMR')

# r-universe (see causal-genomics/mendelian-randomization)
install.packages('TwoSampleMR', repos=c('https://mrcieu.r-universe.dev','https://cloud.r-project.org'))

# Optional: longitudinal / time-varying confounding
install.packages('gfoRmula')
```

Compute requirements:
- Single-mediator + 5000 bootstrap sims: minutes on a laptop.
- HIMA on ~500k CpGs, n=500: 30-60 min with `parallel=TRUE, ncore=8`.
- BAMA Bayesian MCMC: 1-6 hours depending on chain length.
- medDML cross-fitted random forests: 10-30 min for n=2000.

## Quick Start

Tell the AI agent what kind of mediation question is being asked:
- "Test whether expression of GENE_X mediates the effect of rs12345 on disease risk; include sensitivity to unmeasured confounding"
- "Run 4-way decomposition for treatment-mediator-outcome with exposure-mediator interaction"
- "Screen all 450k CpGs for mediators of smoking-lung-cancer association using HIMA2"
- "Run MR-mediation to estimate the fraction of the BMI-CHD effect that goes through LDL"
- "Use double-ML mediation to estimate ACME with random-forest nuisance models"
- "Compute the mediational E-value for an ACME of 0.04 with 95% CI [0.01, 0.07]"

## Example Prompts

### Single-Mediator Observational

> "I have individual-level data on genotype, gene expression, and binary disease status with covariates age/sex/PCs. Test whether expression mediates the genotype-disease association and report ACME with 95% CI and proportion mediated."

> "Run mediation analysis for SNP rs7412 with LDL cholesterol as the mediator and Alzheimer's disease as the outcome; include Imai sensitivity analysis."

### 4-Way Decomposition with Interaction

> "Decompose the smoking-COPD effect into CDE, PIE, INTref, INTmed using FEV1 as the mediator and CMAverse; outcome is binary with rare-disease prevalence around 5%."

> "Fit a regression-based 4-way decomposition with exposure-mediator interaction and exposure-set-to-0 reference."

### High-Dimensional EWAS / Transcriptome

> "Run HIMA2 on the Illumina EPIC methylation matrix to identify CpGs mediating the prenatal-smoking-birthweight relationship; control FDR at 0.05."

> "I have RNA-seq counts for 20k genes, n=200 subjects, and a binary case/control outcome. Use HIMA-binomial to find transcripts that mediate exposure-disease."

> "Run BAMA with Bayesian shrinkage on 5000 candidate protein mediators."

### MR-Mediation

> "Run two-step MR for BMI -> diastolic blood pressure -> coronary artery disease using independent instruments at each step; apply Steiger filter."

> "Use MVMR-mediation to estimate the direct effect of LDL on CHD adjusting for HDL; report conditional F for each exposure and the indirect path through HDL."

### Longitudinal / Time-Varying Confounding

> "Run g-formula mediation via CMAverse for treatment exposure with time-varying confounders measured at three timepoints."

### Sensitivity

> "Compute the mediational E-value for my ACME of 0.05 (95% CI 0.02-0.08) on a binary outcome."

> "Run medsens on this mediate result and tell me whether rho_crit is above 0.3."

## What the Agent Will Do

1. Identify the analytic regime: single vs high-D mediators; observational vs IV-based; continuous / binary / survival outcome; presence of suspected E-M interaction.
2. Choose method: `mediation` for default observational; `CMAverse` for interaction or non-trivial outcome models; `HIMA` / `HIMA2` / `BAMA` for high-D; MR-mediation when IVs are stronger evidence than measured confounder set; `medDML` for doubly-robust under rich confounders.
3. Fit mediator and outcome models with appropriate link functions (linear, logistic, Cox).
4. Bootstrap or analytic CIs (>= 1000 sims exploratory, >= 5000 publication).
5. Report ACME / ADE / total effect / proportion mediated AND a sensitivity result (Imai rho_crit or mediational E-value).
6. For high-D: apply BH FDR <= 0.05 and report screened mediator count.
7. For MR-based: verify conditional F > 10, Steiger filter, and pleiotropy diagnostics.
8. Flag failure modes encountered (exposure-induced confounder, weak instruments, sparse outcomes).

## Tips

- Sequential ignorability is untestable; ALWAYS include sensitivity (Imai rho or E-value). Reviewers reject mediation papers without it.
- HIMA1 (`hima_classic`) screens by beta only and reverses the screening regime relative to HIMA2; HIMA2 (`hima`) screens by alpha-beta -- HIMA2 is more powerful for mediators with strong exposure-effect (alpha) but weak outcome-effect (beta). Pin `HIMA >= 2.3.0` for the formula interface used in this skill's examples; 2.2.x users must use `hima_classic()` style with positional arguments and lose the auto-detected outcome family.
- HIMA outcome family is limited to gaussian and binomial in the default `hima()`; survival needs `hima_cox`, count needs `hima_pois`.
- Bootstrap iterations: 1000 is the floor; 5000 for publication. Doubling sims halves Monte-Carlo CI noise.
- BCa CI caveats: pass `boot.ci.type='bca'` (lowercase) NOT `'BCa'`; BCa requires more sims than percentile (recommend >= 5000) and can fail to converge when the acceleration estimate is unstable near boundary cases. Fall back to percentile CIs in that case and document.
- Proportion-mediated > 0.2 is a soft convention -- absolute ACME size matters more, especially for binary outcomes where pm is unstable near total = 0.
- Difference and product-of-coefficients identical only for linear-Gaussian; for any non-linear outcome use the counterfactual ACME from `mediate()` or CMAverse, not the hand-computed product.
- Binary rare outcomes (<= 10%) allow OR-based 4-way decomposition (Valeri & VanderWeele 2013); common outcomes require RR-scale or marginal-effects scale.
- For exposure-induced confounder of M-Y, switch from natural to interventional indirect effects via `CMAverse::cmest(estimation='msm')` or `gfoRmula`.
- Compare observational ACME against MR-mediation when possible; convergent evidence is high-confidence, divergent evidence flags methodological issues.
- For MVMR-mediation, conditional F > 10 for each exposure is strict; weak instruments inflate direct effect and shrink indirect.
- `medsens()` requires linear or probit outcome models -- refit with `family=binomial(link='probit')` for sensitivity on binary outcomes.
- Cell composition is the canonical confounder in EWAS mediation; include estimated cell proportions (Houseman / RPC) in `COV.XM` and `COV.MY`.
- BAMA when-to-prefer: choose BAMA over HIMA2 when (1) strong prior information on mediator effects is available, (2) the candidate panel is moderately-sized (~5k mediators), and (3) compute budget allows 1-6 h MCMC; otherwise HIMA2 is faster with comparable FDR control.
- Measurement-error correction: when the mediator is measured with error (reliability r < 0.9), the indirect effect is attenuated; apply regression calibration (Carroll 2006 Measurement Error in Nonlinear Models) by replacing the observed mediator with its conditional expectation given the exposure and covariates, OR run a sensitivity analysis at fixed reliability r = 0.7 (Valeri & VanderWeele 2014 Biometrics 70:268).
- Mediational E-value formula details: convert ACME to a risk-ratio bound (`acme_rr = exp(ACME)` on the log scale for continuous, or via VanderWeele's marginal RR for binary), then E = RR + sqrt(RR * (RR - 1)); the corresponding lower-bound RR gives the E-value for the CI. `EValue::evalues.OLS()` automates this for linear outcomes.

## Cohort Considerations

Large biobank cohorts have well-characterized strengths and weaknesses for mediation analysis:
- UK Biobank: largest n, rich phenotyping, individual-level access via application. EWAS mediation feasible at ~50k subset with EPIC arrays. Survival outcomes via NHS linkage.
- FinnGen: well-suited to MR-mediation (summary-stat releases) and time-to-event outcomes (universal-registry follow-up); individual-level mediator measurements limited.
- ALSPAC (Avon): longitudinal mediator measurements ideal for time-varying mediation (g-formula, MSM); smaller n caps high-D screening power.

## Related Skills

- causal-genomics/mendelian-randomization - IV-based causal inference; foundation for MR-mediation
- causal-genomics/colocalization-analysis - Confirm shared causal variant before causal mediation
- causal-genomics/fine-mapping - Identify the causal variant driving the exposure
- methylation-analysis/differential-cpg-testing - Per-CpG inputs for HIMA EWAS mediation
- differential-expression/deseq2-basics - Expression inputs for eQTL mediation
- multi-omics-integration/mofa-integration - Multi-layer mediator construction
- population-genetics/association-testing - GWAS summary statistics for MR-mediation
- clinical-biostatistics/effect-measures - Risk-ratio / odds-ratio scales for binary outcomes
- machine-learning/model-validation - Cross-fitting and sample splitting for medDML

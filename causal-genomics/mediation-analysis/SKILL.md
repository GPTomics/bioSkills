---
name: bio-causal-genomics-mediation-analysis
description: Decompose total effects into direct and indirect paths through mediators using mediation, CMAverse 4-way, HIMA/HIMA2 high-dimensional, BAMA, two-step / MVMR mediation, or double-ML medDML. Use when testing whether a molecular phenotype (expression, methylation, protein) mediates a treatment-outcome relationship, decomposing exposure-mediator interaction via VanderWeele 4-way, screening high-dimensional EWAS mediators, or running MR-based mediation when sequential ignorability is implausible.
tool_type: mixed
primary_tool: mediation
---

## Version Compatibility

Reference examples tested with: R 4.3+, mediation 4.5.0+, CMAverse 0.1.0+ (GitHub `BS1125/CMAverse`), HIMA >= 2.3.0 (CRAN), bama 1.3+, causalweight 1.0.5+ (medDML), MVMR 0.4+, TwoSampleMR 0.6+, EValue 4.1+, gesttools 1.3+, ipw 1.0.11+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- HIMA must be pinned at `>= 2.3.0` for the code patterns below; the formula interface `hima(formula, data.pheno, data.M, mediator.type, penalty, ...)` was introduced in 2.3.0. Users on HIMA 2.2.x must instead call `hima_classic()` with positional arguments (`X=`, `Y=`, `M=`, `COV.XM=`, `Y.family=`) and lose the formula interface; HIMA 2.2.x is NOT API-compatible with the examples here.
- `hima_classic()` is still available in 2.3+ for the original Zhang 2016 SIS+MCP pipeline; use it only to reproduce 2016-2021 papers.

If code throws an error, introspect the installed package (`?hima`, `args(cmest)`) and adapt the example to match the actual API rather than retrying.

# Mediation Analysis

**"Does expression of GENE_X mediate the SNP-to-disease effect?"** -> Decompose the total effect of a treatment (genotype, exposure) on an outcome into direct and indirect paths through one or more mediators, with explicit handling of exposure-mediator interaction, sensitivity to unmeasured confounding, and high-dimensional mediator screening.

- R (single-mediator, observational, sequential-ignorability assumed): `mediation::mediate(med_model, out_model, treat='X', mediator='M', boot=TRUE, sims=5000)`
- R (4-way decomposition with exposure-mediator interaction): `CMAverse::cmest(...EMint=TRUE, estimation='paramfunc', inference='bootstrap', nboot=1000)`
- R (high-dimensional / EWAS mediators): `HIMA::hima(Y ~ X + covariates, data.pheno, data.M, mediator.type='gaussian', penalty='DBlasso')` (modern v2.3+ formula interface)
- R (MR-based mediation): two-step `TwoSampleMR` with independent instruments OR `MVMR::ivw_mvmr` for joint direct effect
- R (doubly-robust double-ML): `causalweight::medDML(y, d, m, x)`

Sequential ignorability (no unmeasured confounder of treatment-mediator, mediator-outcome, treatment-outcome) is the single load-bearing assumption of observational mediation and is fundamentally untestable. Every report should include a sensitivity result (Imai's rho via `medsens()` or a mediational E-value).

## Algorithmic Taxonomy

| Method | Framework | Handles E-M interaction | High-D mediators | Min n | Fails when |
|--------|-----------|--------------------------|------------------|-------|------------|
| Baron-Kenny (1986) | Additive regression-based product/difference | No | No | ~100 | Any non-linearity, interaction, or binary outcome; deprecated for causal inference |
| Imai mediation R (Imai 2010 Psychol Methods) | Counterfactual ACME/ADE with bootstrap | Yes (via interaction term in outcome model) | No | ~200 | Sequential ignorability violated; exposure-induced M-Y confounder; rare binary outcome with logistic outcome model |
| VanderWeele 4-way (VanderWeele 2014 Epidemiology 25:749; 2015 OUP) | CDE + PIE + INTref + INTmed decomposition | Native | No | ~300 | Without interaction term reduces to standard mediation; binary outcome needs rare-disease assumption |
| CMAverse (Shi 2021 Epidemiology 32:e20) | 6 estimators: regression (`rb`), weighting (`wb`), IORW (`iorw`), natural effect models (`ne`), MSM (`msm`), g-formula (`gformula`) | Yes | No (single M, or M-vector) | ~300 | Estimator-specific; `wb` fails with rare exposure; `msm` needs censoring weights for survival |
| HIMA1 / hima_classic (Zhang 2016 Bioinformatics 32:3150) | SIS screen by beta (M->Y) + MCP penalty | No | Yes (up to ~10k) | ~150 + p>>n | Misses mediators with strong alpha and weak beta; screening-step false-negatives |
| HIMA2 / hima (Perera 2022 BMC Bioinformatics 23:296) | SIS screen by alpha*beta (indirect effect) + MCP | No (linear by default) | Yes | ~150 | Outcome family limited (gaussian/binomial); HIMA-Cox for survival; HIMA-Pois for count |
| HILAMA (Zhang 2025) | High-D mediation with latent confounders | No | Yes (>= 100k) | ~500 | Newer; benchmarks evolving; requires latent-factor specification |
| BAMA (Song 2020 Biostatistics) | Bayesian high-D continuous shrinkage | No | Yes (~5k) | ~200 | Slow MCMC; prior sensitivity for very weak mediators |
| Two-step MR / network MR (Burgess 2017 Eur J Epidemiol 32:377) | IV-based at each step with INDEPENDENT instruments | Implicit (no interaction modeling) | One mediator at a time | Large summary-stat samples | Same SNP used for E and M (violates exclusion); horizontal pleiotropy; Steiger reversal of M->E direction |
| MVMR-mediation (Carter & Sanderson 2021 Eur J Epidemiol 36:465-478) | Total minus direct via MVMR | Implicit | Single mediator | Large GWAS samples for both E and M | Conditional F < 10 for either exposure; correlated instruments |
| medDML (Farbmacher 2022 Econometrics J 25:277) | Double-debiased ML, doubly-robust | Limited (depends on learner) | Moderate (sparsity-friendly) | ~500 | Severe overlap violations; cross-fitting variance with small n |

Methodology evolves; verify against the current CMAverse vignette and the Steen / Vansteelandt natural-effects-model literature before locking analytic choices. Difference-in-coefficients and product-of-coefficients give identical estimates in fully linear-Gaussian models but DIVERGE for any non-linear outcome model (logistic, Cox, Poisson); the counterfactual ACME from `mediation::mediate()` is the correct quantity for non-linear outcomes.

## 4-Way Decomposition Framework

VanderWeele 4-way explicitly separates effects from exposure-mediator interaction:

```
Total Effect = CDE + INTref + INTmed + PIE
```

| Component | Meaning | Active when |
|-----------|---------|-------------|
| CDE | Controlled direct effect (with mediator fixed at reference level) | E directly affects Y |
| INTref | Interaction-reference -- needs interaction AND exposure | E*M interaction with mediator at reference |
| INTmed | Mediated interaction -- needs interaction AND exposure AND mediation | E shifts M which then interacts with E |
| PIE | Pure indirect effect (older "mediation" quantity) | E shifts M which shifts Y additively |

Without an exposure-mediator interaction term, INTref = INTmed = 0 and the decomposition collapses to CDE + PIE (= ADE + ACME). With interaction present, traditional ACME mixes PIE and INTmed; the 4-way separation is the only framework that disentangles them. Most epidemiology applications include the interaction term and report all four components (Valeri & VanderWeele 2013 Psychol Methods).

## Decision Tree by Scenario

| Scenario | Recommended pipeline |
|----------|---------------------|
| Observational, single measured mediator, no plausible E-M interaction, continuous outcome | `mediation::mediate()` with `boot=TRUE, sims=5000`; always run `medsens()` |
| Observational, single mediator, suspected E-M interaction, any outcome family | `CMAverse::cmest(..., EMint=TRUE)` -> read CDE, PIE, INTref, INTmed |
| Observational, BINARY outcome, rare disease (< 10%) | `cmest(yreg='logistic', EMint=TRUE, casecontrol=FALSE)` -- OR-based 4-way decomposition is valid under rare-disease |
| Observational, survival outcome | `cmest(yreg='coxph')` OR `HIMA::hima_cox` for high-D; report HRs |
| High-D mediators (EWAS, transcriptome-wide), continuous outcome | `HIMA::hima(formula, data.pheno, data.M, mediator.type='gaussian', penalty='DBlasso')`; report `sigcut` (FDR threshold, default 0.05) |
| High-D mediators with latent confounding (very-high-D EWAS) | `HILAMA` (2025) |
| High-D mediators with Bayesian shrinkage (small n, ~5k features) | `bama::bama()` |
| Strong genetic IVs for exposure available, single mediator with own IVs | Two-step MR with independent instruments + Steiger filter on mediator |
| Both E and M have IVs but instruments are weak / correlated | MVMR-mediation with conditional F > 10 each |
| Observational with rich confounder set, want doubly-robust estimate | `causalweight::medDML` (double-debiased ML) |
| Longitudinal with time-varying confounding | g-formula via `CMAverse::cmest(estimation='gformula')` OR `gfoRmula` package |
| Exposure-induced confounder of M-Y exists | Interventional indirect effects (Vansteelandt & Daniel 2017); `CMAverse::cmest(estimation='msm')` |

## Sequential Ignorability and Why It Always Needs Sensitivity

Observational mediation requires three no-unmeasured-confounding assumptions. The third (M-Y unmeasured confounder, after conditioning on E) is the most common violator in genomic mediation because biological confounders (cell composition, batch effects, technical mediators) frequently affect both M and Y.

### Sequential ignorability untestable

**Trigger:** Always, by design.

**Mechanism:** No statistical test can detect an unmeasured confounder of M-Y. Bootstrap CIs assume the assumption holds; they do NOT propagate uncertainty about it.

**Symptom:** Significant ACME with no sensitivity reported -> reviewer rejects.

**Fix:** Report at least one of:
- Imai's rho sensitivity: `medsens(med_result, rho.by=0.05, sims=1000)`; the critical rho where ACME crosses 0; |rho_crit| > 0.3 is "reasonably robust" (Imai 2010), |rho_crit| < 0.1 is highly sensitive.
- Mediational E-value (Smith & VanderWeele 2019 Epidemiology 30:835): minimum risk-ratio strength of an unmeasured confounder needed to nullify the observed indirect effect; computed via `EValue::evalues.OLS()` for linear outcomes or by-hand from ACME risk ratio bounds.
- Reporting BOTH rho-based and E-value sensitivity is standard for high-stakes claims.

### Methods-Section Defense of Sequential Ignorability

Template sentence for the methods write-up: "We assumed sequential ignorability conditional on {age, sex, ancestry PCs, cell composition, batch, smoking}. Robustness was assessed via Imai rho_crit at the ACME contrast (`medsens`, sims = 1000) and the mediational E-value on the risk-ratio scale (Smith & VanderWeele 2019 Epidemiology 30:835)."

Quantitative interpretation thresholds:
- rho_crit: |rho_crit| > 0.3 robust; 0.1-0.3 moderately sensitive; < 0.1 highly sensitive (Imai 2010).
- E-value: E > 2 robust to plausible biological confounding; 1.5-2 moderate; < 1.5 fragile (Smith & VanderWeele 2019).

For high-stakes claims (clinical, drug-target, regulatory submissions) report BOTH rho_crit and the mediational E-value; for exploratory work either alone suffices.

### Exposure-induced M-Y confounder

**Trigger:** A covariate L sits between E and Y, AND is affected by E, AND confounds M-Y.

**Mechanism:** Standard regression-based mediation cannot adjust for L without blocking part of the indirect effect (collider stratification bias). Adjusting biases CDE; not adjusting biases ACME.

**Symptom:** Sensitivity to confounder set; ACME flips sign when L is added vs removed.

**Operational identification:** From the DAG, L is a covariate of M and Y that is also affected by E. Vansteelandt 2009 (Stat Sci 24:471) gives the criterion: if L is adjusted, part of the indirect E -> L -> M -> Y pathway is blocked; if L is not adjusted, L confounds the M-Y leg. Both are wrong under natural-effects; the natural indirect effect is simply not identified.

**Fix:** Switch to **interventional indirect effects** (Vansteelandt & Daniel 2017 Epidemiology 28:258), NOT natural indirect effects. Use `CMAverse::cmest(estimation='msm')` with stabilized inverse-probability weights (yields the randomized-interventional analogue), `gfoRmula` (parametric g-formula), or randomized indirect effects (Lin 2017 Biometrics 73:1109). The interventional indirect is identified under weaker assumptions than the natural indirect.

### HIMA covariate or data.pheno error

**Trigger:** `data.pheno` contains factor columns with NA, or formula references columns missing from `data.pheno`.

**Mechanism:** HIMA v2.3+ uses a formula interface and constructs the design matrix internally from `data.pheno`; missing values or unparseable formulas surface as cryptic `glmnet` errors.

**Symptom:** Pipeline fails inside `hima()` with a non-obvious `storage.mode` or `model.matrix` error.

**Fix:** Pre-clean `data.pheno` (drop NA rows for the variables in the formula; convert factors with `factor()`; ensure all RHS variables exist as columns). Example:
```r
dat <- na.omit(dat[, c('outcome', 'exposure', 'age', 'sex', 'batch', 'pc1', 'pc2')])
dat$batch <- factor(dat$batch)
result <- hima(outcome ~ exposure + age + sex + batch + pc1 + pc2,
               data.pheno=dat, data.M=M_matrix, mediator.type='gaussian')
```

### HIMA mediator-type vs outcome-type mismatch

**Trigger:** Survival outcome (`Surv()` on LHS) with `mediator.type='compositional'` chosen against text mediator panel; or count mediators passed as `'gaussian'`.

**Mechanism:** HIMA v2.3+ auto-detects outcome family from the LHS of `formula` (continuous, binary, survival, count); `mediator.type` is set for the mediator data only (`'gaussian'`, `'negbin'`, `'compositional'`). Mismatching mediator-type to the actual mediator distribution biases the screening step.

**Symptom:** Hazard / rate ratios for indirect effects look implausible; many "significant" mediators fail replication.

**Fix:** Set `mediator.type='gaussian'` for continuous (e.g., methylation beta, log-CPM expression), `'negbin'` for raw count (RNA-seq), `'compositional'` for relative-abundance microbiome. Verify by `?hima` in the installed version since the catalogue of mediator types has expanded across releases.

### Bootstrap iterations too low

**Trigger:** `sims=100` or `sims=500` in early exploration left in for the final report.

**Mechanism:** ACME CIs from bootstrap have Monte-Carlo error that scales as 1/sqrt(sims); at sims=500 the 95% CI bounds have ~5% MC noise, enough to flip the conclusion at the boundary.

**Symptom:** Re-running `mediate()` with a different `set.seed()` gives substantially different CI bounds.

**Fix:** `sims=1000` minimum for any reported result; `sims=5000` for publication; `sims=10000` if proximity to zero matters. BCa CIs (`boot.ci.type='bca'` in `mediate()`) are slightly more accurate than percentile CIs near zero but require more sims for stability.

### Two-step MR instrument independence

**Trigger:** Same set of SNPs used as instruments for E in step 1 and for M in step 2.

**Mechanism:** If a SNP affects both E and M, the M-instrument violates exclusion restriction (the SNP-Y association is not exclusively through M). Estimates are biased toward the direct effect.

**Symptom:** Two-step MR shows large indirect effect; replacing M-instruments with non-overlapping SNPs makes it vanish.

**Fix:** Apply Steiger filter on the mediator: keep only SNPs where the SNP-M F-statistic exceeds SNP-E F-statistic (or where SNP explains more variance in M than E). For MVMR-mediation, require conditional F > 10 for each exposure independently (Sanderson 2019 IJE 48:713).

### Difference vs product of coefficients diverge for non-linear outcomes

**Trigger:** Binary or survival outcome modeled with logistic / Cox.

**Mechanism:** Difference = total - direct; product = alpha * beta. Equivalent under linear-Gaussian; diverge under any link function. Counterfactual ACME from `mediation::mediate()` is the correct quantity; hand-computed product-of-coefficients on logistic output is biased except under rare-disease.

**Fix:** Report only counterfactual ACME (Imai or CMAverse). For OR-based 4-way decomposition on rare outcomes (<= 10%), Valeri & VanderWeele 2013 formulas apply; for common outcomes use risk-ratio scale or marginal effects rather than ORs.

### Required Reporting for Publication

| Component | Required |
|-----------|----------|
| ACME estimate + 95% CI (BCa preferred) | Yes |
| ADE + 95% CI | Yes |
| Total effect | Yes |
| Proportion mediated | Yes when total > effect-size threshold |
| Bootstrap method + sims | percentile / BCa; min 1000, recommend 5000 |
| Sequential ignorability sensitivity | rho_crit (medsens) OR mediational E-value |
| Exposure-mediator interaction test | Coefficient + p; 4-way decomposition if significant |
| Confounder set justification | DAG description |
| Mediator measurement reliability | Cite |
| Sample size + missing-data handling | Yes |
| Mediator / exposure scale | Standardized? log? raw? |

Reference: AGReMA-Mediation guideline (Lee 2021 BMJ 372:n122) and MacKinnon 2008 Introduction to Statistical Mediation Analysis.

## Reconciliation: Observational vs MR Mediation

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| Observational ACME significant; MR-mediation null | Unmeasured M-Y confounding inflated observational estimate; OR weak IVs in MR | Re-run observational with `medsens()`; if rho_crit < 0.1, trust MR null |
| Observational ACME null; MR-mediation significant | Measurement error in M attenuated observational estimate | Trust MR (regression-dilution-free) IF instruments pass Steiger and pleiotropy tests (MR-Egger intercept, MR-PRESSO) |
| Both significant with same sign | Convergent evidence | High-confidence mediation; report effect size from the more-conservative estimate |
| Both significant with opposite signs | At least one is biased; revisit confounder structure and IV assumptions | Do not pool; investigate via cross-method sensitivity |

**Operational rule for high-stakes claims (clinical / drug-target mediation):** Require (1) significant observational ACME, (2) Imai rho_crit > 0.2 OR mediational E-value > 1.5, (3) directionally consistent MR-mediation result OR documented absence of valid instruments. Single-method mediation claims should be reported as exploratory.

## Anticipated Reviewer Pushback

| Pushback | Standard response |
|----------|-------------------|
| "Sequential ignorability?" | Imai rho_crit reported via `medsens`; mediational E-value reported on the risk-ratio scale |
| "Exposure-induced confounder of M-Y?" | DAG drawn; if L present, switch to `CMAverse::cmest(estimation='msm')` for interventional indirect effect (Vansteelandt & Daniel 2017) |
| "Why this bootstrap method?" | BCa with sims=5000 for publication; percentile fallback when BCa fails to converge (acceleration estimate unstable at boundary) |
| "Why didn't you do MR-mediation?" | If valid IVs for E and M exist: two-step MR or MVMR-mediation done (see code below); if not, documented absence of trans-instruments |
| "Mediator measured with error?" | Regression calibration (Carroll 2006 Measurement Error in Nonlinear Models) OR sensitivity analysis assuming reliability r = 0.7 (Valeri & VanderWeele 2014) |
| "Why HIMA2 not BAMA?" | HIMA2 = frequentist + FDR control + faster; BAMA = Bayesian when prior information is available; sample-size justification given against simulation rule-of-thumb |
| "Proportion mediated unstable?" | When |total| < 2*SE(total), proportion-mediated CI is unreliable (denominator near zero); report indirect effect alone with absolute effect size |

## Quantitative Thresholds

| Quantity | Threshold | Source / Rationale |
|----------|-----------|--------------------|
| `sims` (bootstrap iterations) | >= 1000 exploratory, >= 5000 publication | Imai 2010; MC error scales 1/sqrt(sims) |
| Proportion mediated -- meaningful | > 0.2 | Convention; weak guideline only -- effect size in absolute terms matters more (MacKinnon 2008) |
| Proportion mediated -- "most of the effect" | > 0.5-0.8 | Convention |
| Imai rho_crit -- robust | > 0.3 | Imai 2010 Psychol Methods 15:309 |
| Imai rho_crit -- sensitive | < 0.1 | Same |
| Mediational E-value -- robust | > 2.0 (working convention; the original Smith & VanderWeele 2019 E-value framework does not prescribe a specific cutoff -- magnitude is context-dependent) | Smith & VanderWeele 2019 Epidemiology 30:835 |
| HIMA FDR cutoff | BH FDR < 0.05 | Default; report q-values not raw p |
| MVMR conditional F per exposure | > 10 each | Sanderson 2019 IJE 48:713 |
| Two-step MR -- F for both stages | > 10 each | Burgess weak-instrument convention |
| Sample size -- single-mediator (Imai) | >= 200 for stable bootstrap | Simulation rule-of-thumb |
| Sample size -- HIMA EWAS | >= 150 with p_mediators up to ~10k | Zhang 2016 simulations |
| Rare-outcome cutoff for OR-based 4-way | outcome prevalence <= 10% | Valeri & VanderWeele 2013 |

## Working Code Patterns

### Single-Mediator Observational with Sensitivity

**Goal:** Decompose a genotype-disease effect via measured molecular mediator with explicit sensitivity to unmeasured M-Y confounding.

**Approach:** Fit mediator and outcome models, bootstrap ACME/ADE, run `medsens()` for Imai rho-based sensitivity.

```r
library(mediation)

med_model <- lm(expression ~ genotype + age + sex + pc1 + pc2 + pc3, data=dat)
out_model <- glm(disease ~ genotype + expression + age + sex + pc1 + pc2 + pc3,
                 data=dat, family=binomial)

med_result <- mediate(med_model, out_model,
                      treat='genotype', mediator='expression',
                      boot=TRUE, sims=5000, boot.ci.type='bca')
summary(med_result)

sens <- medsens(med_result, rho.by=0.05, effect.type='indirect', sims=1000)
summary(sens)
```

`d0`, `z0`, `n0`, `tau.coef` slots return ACME, ADE, proportion mediated, total effect.

### 4-Way Decomposition with Exposure-Mediator Interaction

**Goal:** Separate CDE, PIE, INTref, INTmed when exposure-mediator interaction is biologically plausible (e.g., gene-environment interaction modifying mediator effect).

**Approach:** Use CMAverse regression-based estimator with `EMint=TRUE`; bootstrap CIs.

```r
library(CMAverse)

result_4way <- cmest(
  data=dat, model='rb',
  outcome='disease', exposure='genotype', mediator='expression',
  basec=c('age','sex','pc1','pc2'),
  EMint=TRUE,
  mreg=list('linear'), yreg='logistic',
  astar=0, a=1, mval=list(0),
  estimation='paramfunc', inference='bootstrap', nboot=1000
)
summary(result_4way)
```

CMAverse reports the 4-way decomposition (Vanderweele 2014): for continuous outcomes the components are `cde`, `intref`, `intmed`, `pnie` (or `pie`), `te`, `pm`; for non-continuous outcomes (logistic / Cox / Poisson) the ratio versions `Rcde`, `Rpnde`, `Rtnde`, `Rpnie`, `Rtnie` are reported. When `EMint=TRUE`, additional proportion-attributable-to-interaction terms (`int`, `pe`) are included. Verify column names with `summary(result)$results` in the installed CMAverse version, since naming has evolved.

### High-Dimensional EWAS Mediation (HIMA2)

**Goal:** Among thousands of candidate CpG mediators, identify those mediating an exposure-outcome effect with FDR control.

**Approach:** HIMA v2.3+ uses a formula interface and auto-detects outcome family (Gaussian / binomial / Cox / Poisson). Screening + MCP/DBlasso penalisation + joint significance with BH; the `sigcut` argument controls the FDR threshold (default 0.05).

```r
library(HIMA)

dat <- na.omit(dat[, c('outcome', 'exposure', 'age', 'sex', 'cell_pc1', 'cell_pc2')])
M_matrix <- as.matrix(beta_values)

result <- hima(
  outcome ~ exposure + age + sex + cell_pc1 + cell_pc2,
  data.pheno=dat,
  data.M=M_matrix,
  mediator.type='gaussian',          # 'negbin' for count, 'compositional' for microbiome
  penalty='DBlasso',                  # default; alternatives 'MCP', 'SCAD', 'lasso'
  scale=TRUE,
  sigcut=0.05,
  parallel=TRUE, ncore=8, verbose=TRUE
)
# result is a data.frame of significant mediators below sigcut
```

For survival outcomes wrap the LHS as `Surv(time, status)`; HIMA auto-routes to Cox. The old `hima_classic()` (Zhang 2016 original) is still exported but screens by beta only and misses mediators with strong alpha + weak beta -- prefer the wrapper `hima()` unless reproducing a 2016-2021 paper.

For highly-correlated mediators (CpG-island clusters, gene-module co-expression): HIMA uses joint significance with BH-FDR on max(p_alpha, p_beta) and handles correlation only weakly. Within `hima()`, set `penalty='MCP'` for stronger correlation handling; alternatively pre-reduce the mediator panel by principal components or by clustering correlated mediators and screening the cluster centroid (VanderWeele & Vansteelandt 2014 Epidemiol Methods 2:95).

### Time-Varying Mediation Methods

When the mediator is measured at multiple timepoints (or exposure varies over time), natural-effects estimands are not identified; switch to one of:

- g-formula (parametric or Monte Carlo): `gfoRmula::gformula_continuous_eof()`; `CMAverse::cmest(estimation='gformula')`
- g-estimation of a structural nested mean model: `gesttools::gestSingle()` / `gestMultiple()`
- Marginal structural model with stabilized IPTW: `ipw::ipwtm()` followed by `glm(..., weights=sw)`
- Sequential mediation for K timepoints: VanderWeele & Tchetgen Tchetgen 2017 JRSSB 79:917

Choose based on the experimental structure:
- >= 3 timepoints required for g-methods to identify time-varying indirect effects
- Longitudinal mediator measurement at EACH timepoint is required (not just baseline)
- MSM is preferred when treatment is binary and time-varying; g-formula when continuous
- Sequential mediation when the causal ordering of multiple mediators is known and stable across time

### MR-Mediation: Two-Step vs MVMR-Mediation

Decision tree:
- Independent instrument sets available for E and M -> two-step MR (Burgess 2017 Eur J Epidemiol 32:377)
- E and M share instruments (common in cis-eQTL / cis-pQTL mediator cases) -> MVMR-mediation (Carter & Sanderson 2021 Eur J Epidemiol 36:465)
- Both feasible -> report both (triangulation)

Two-step code sketch:

```r
library(TwoSampleMR)
exp_E <- extract_instruments('ieu-a-2', clump=TRUE)
m_E <- extract_outcome_data(exp_E$SNP, 'ieu-b-30')
dat_EM <- harmonise_data(exp_E, m_E)
mr_EM <- mr(dat_EM)

exp_M <- extract_instruments('ieu-b-30', clump=TRUE)
exp_M_indep <- exp_M[!exp_M$SNP %in% exp_E$SNP, ]
exp_M_indep <- steiger_filtering(exp_M_indep)
out_M <- extract_outcome_data(exp_M_indep$SNP, 'ieu-a-7')
dat_MY <- harmonise_data(exp_M_indep, out_M)
mr_MY <- mr(dat_MY)
```

Indirect effect = beta_EM * beta_MY (product of coefficients). CI via delta method or parametric bootstrap of the joint (beta_EM, beta_MY) distribution. Steiger filter on M-instruments is mandatory to ensure the M -> Y direction (not Y -> M).

### MR-Mediation: Total Minus Direct via MVMR

**Goal:** Estimate the proportion of a genetic-instrument-identified causal effect that flows through a mediator, using independent IVs for E and (E + M).

**Approach:** Univariable MR for total E->Y; MVMR for direct E->Y conditional on M; indirect = total - direct via delta-method CI.

```r
library(TwoSampleMR); library(MVMR)

total <- mr_ivw(harmonised_data_E_Y)
mvmr_dat <- format_mvmr(BXGs=cbind(beta_E, beta_M),
                        BYG=beta_Y, seBXGs=cbind(se_E, se_M), seBYG=se_Y, RSID=snps)
fstat <- strength_mvmr(mvmr_dat, gencov=0)
mvmr_fit <- ivw_mvmr(mvmr_dat)
direct <- mvmr_fit$coef[1, 'Estimate']
direct_se <- mvmr_fit$coef[1, 'Std. Error']

indirect <- total$b - direct
indirect_se <- sqrt(total$se^2 + direct_se^2)
indirect_ci <- indirect + c(-1.96, 1.96) * indirect_se
```

Require `fstat` conditional F > 10 for both E and M independently. If `fstat < 10`, use Q-statistic-adjusted IVW (`qhet_mvmr`) or report the result as weak-instrument-limited.

### Double-ML Doubly-Robust Mediation

**Goal:** Avoid model misspecification of both mediator and outcome models via cross-fitted ML nuisance estimators.

**Approach:** `causalweight::medDML` uses random forests (or other learners) with sample splitting to estimate nuisance parameters; final estimator is doubly robust.

```r
library(causalweight)

result_dml <- medDML(
  y=dat$outcome, d=dat$treatment, m=dat$mediator,
  x=as.matrix(dat[, covariates]),
  trim=0.05, order=1
)
```

Reports direct, indirect (via mediator), and total effects with influence-function-based standard errors. Robust to non-linearity and interactions; assumes sequential ignorability still.

### Mediational E-Value for Sensitivity

**Goal:** Report the minimum strength of an unmeasured M-Y confounder required to nullify the observed indirect effect.

**Approach:** Convert ACME and its CI to a risk-ratio scale, then apply VanderWeele E-value formula.

```r
library(EValue)

acme_rr <- exp(med_result$d0)
acme_lower_rr <- exp(med_result$d0.ci[1])
evalues.RR(acme_rr, lo=acme_lower_rr, hi=NULL)
```

For binary outcomes, convert ACME on probability scale to RR; for continuous, use `evalues.OLS()` with the standardized indirect effect. E-value > 2 indicates a confounder would need >2-fold associations with both M and Y to nullify the indirect effect (Smith & VanderWeele 2019).

## Tool Install Notes

| Package | Source | Notes |
|---------|--------|-------|
| mediation | CRAN | `install.packages('mediation')`; actively maintained (Imai group) |
| CMAverse | GitHub | `remotes::install_github('BS1125/CMAverse')`; NOT on CRAN; 6 estimators in one interface |
| HIMA | CRAN | `install.packages('HIMA')`; v2.x renamed `hima()` to HIMA2 -- verify with `?hima` |
| bama | CRAN | `install.packages('bama')`; Bayesian; slow MCMC |
| causalweight | CRAN | `install.packages('causalweight')`; medDML for double-ML mediation |
| EValue | CRAN | `install.packages('EValue')`; for mediational E-values |
| TwoSampleMR | r-universe | See causal-genomics/mendelian-randomization for setup |
| MVMR | r-universe | `remotes::install_github('WSpiller/MVMR')`; for MVMR-mediation |
| gfoRmula | CRAN | For longitudinal / time-varying confounders |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `Error in storage.mode(x) <- "double"` inside `hima()` | NA in `data.pheno` columns referenced by formula, or unconverted factors | `na.omit(data.pheno)` first; ensure all RHS vars in formula are numeric or factor |
| ACME significant, ADE significant, total NOT significant | Suppression / inconsistent mediation | Report transparently; effect partitioning can exceed total in suppression |
| `medsens()` errors on glm outcome | `medsens` requires linear OR probit (not logit) outcome | Refit outcome as `glm(..., family=binomial(link='probit'))` |
| `mediate()` runs forever with binary outcome | `sims=5000` with bootstrap and small n | Use `sims=1000` exploratory; verify model converges first; consider parallel via `parallel='multicore'` |
| CMAverse `cmest()` reports NaN for `pm` | Total effect crosses zero -> proportion ill-defined | Report ACME and TE separately; pm is unstable when |TE| is small |
| Different ACME between `mediation` and CMAverse `rb` | Default `astar/a` levels differ; binary mediator handled differently | Set `astar=0, a=1` explicitly; for binary mediator pass `mval=list(0)` |
| HIMA returns zero significant mediators | Screening too aggressive; or no true mediators | Try `topN=2*sqrt(n)` instead of default; verify with permutation null |
| Two-step MR shows indirect > total | Steiger reversal: M actually causes E; or pleiotropic SNPs | Run MR-Steiger filter; use MR-PRESSO for pleiotropy |
| `medDML` trim removes most data | Severe positivity violation -- few units with overlapping treatment/mediator distributions | Tighten covariate set; check propensity score distributions |

## References

- Baron RM, Kenny DA 1986 J Pers Soc Psychol 51:1173 (original product-of-coefficients)
- Imai K, Keele L, Tingley D 2010 Psychol Methods 15:309 (counterfactual mediation, sequential ignorability)
- VanderWeele TJ 2014 Epidemiology 25:749 (4-way decomposition)
- VanderWeele TJ 2015 Explanation in Causal Inference (OUP) -- canonical textbook
- Valeri L, VanderWeele TJ 2013 Psychol Methods 18:137 (binary outcomes; rare-disease 4-way)
- Vansteelandt S, Daniel RM 2017 Epidemiology 28:258 (interventional / randomized indirect effects)
- Shi B et al 2021 Epidemiology 32:e20 (CMAverse package; 6 estimators)
- Zhang H et al 2016 Bioinformatics 32:3150 (HIMA original)
- Perera C et al 2022 BMC Bioinformatics 23:296 (HIMA2 alpha-beta screening)
- Song Y et al 2020 Biostatistics 21:806 (BAMA Bayesian high-D mediation)
- Burgess S et al 2017 Eur J Epidemiol 32:377 (two-step / network MR)
- Sanderson E et al 2019 IJE 48:713 (MVMR conditional F-statistic)
- Sanderson E et al 2021 IJE 50:1651 (MVMR-mediation tutorial)
- Farbmacher H et al 2022 Econometrics J 25:277 (medDML / double-ML mediation)
- Smith LH, VanderWeele TJ 2019 Epidemiology 30:835 (mediational E-value)

## Related Skills

- causal-genomics/mendelian-randomization - IV-based causal inference; foundation for MR-mediation
- causal-genomics/pleiotropy-detection - MR-mediation instrument validity hinges on pleiotropy diagnostics (MR-Egger intercept, MR-PRESSO)
- causal-genomics/colocalization-analysis - Confirm shared causal variant before causal mediation
- causal-genomics/fine-mapping - Identify the causal variant driving the exposure
- methylation-analysis/differential-cpg-testing - Per-CpG inputs for HIMA EWAS mediation
- differential-expression/deseq2-basics - Expression inputs for eQTL mediation
- multi-omics-integration/mofa-integration - Multi-layer mediator construction
- population-genetics/association-testing - GWAS summary statistics for MR-mediation
- clinical-biostatistics/effect-measures - Risk-ratio / odds-ratio scales for binary outcomes
- machine-learning/model-validation - Cross-fitting and sample splitting for medDML

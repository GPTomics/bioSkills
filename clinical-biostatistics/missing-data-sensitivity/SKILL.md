---
name: bio-clinical-biostatistics-missing-data
description: Implements missing-data sensitivity analyses for confirmatory clinical trials including MMRM under MAR (with Kenward-Roger correction), reference-based multiple imputation (J2R, CR, CIR, LMCF per Carpenter-Roger 2013), Permutt delta-adjustment / tipping-point analysis, pattern-mixture identifying restrictions (CCMV, NCMV, ACMV), and the Cro vs Bartlett variance debate. Use when handling missing primary or secondary endpoint data in regulatory submissions following NRC 2010 and ICH E9(R1).
tool_type: mixed
primary_tool: rbmi
---

## Version Compatibility

Reference examples tested with: R `mmrm` 0.3+ (Roche/openpharma), R `rbmi` 1.5+ (Roche/Bayer via insightsengineering), R `mice` 3.16+, R `mitools` 2.4+, Python `sklearn` 1.4+, `statsmodels` 0.14+.

Before using code patterns, verify installed versions match. If versions differ:
- R: `packageVersion('<pkg>')` then `?function_name`
- Python: `pip show <package>` then `help(module.function)`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Missing Data Sensitivity for Confirmatory Trials

**"Handle missing data in a confirmatory clinical trial"** -> Pre-specify the missing-data assumption per ICH E9(R1); execute the primary analysis under the chosen assumption (typically MAR via MMRM or MI); run clinically-articulable MNAR sensitivity analyses (reference-based MI per Carpenter-Roger 2013); report the tipping delta that would overturn the conclusion (Permutt 2016).

## The Foundation -- NRC 2010 and ICH E9(R1)

**The U.S. National Research Council Panel** ("The Prevention and Treatment of Missing Data in Clinical Trials," 2010; chaired by Roderick Little; Little, D'Agostino, Cohen et al 2012 *NEJM* 367:1355): 18 recommendations grouped as prevention (Recs 1-7), analysis (Recs 8-14), sensitivity (Recs 15-18).

Key recommendations:

- **Rec 10/11:** explicitly REJECT LOCF and BOCF as default; they are biased even under MCAR
- **Rec 13:** endorses WGEE (weighted GEE) for marginal estimands
- **Rec 15:** "examining sensitivity to assumptions about the missing-data mechanism should be a mandatory component of reporting"

**ICH E9(R1) (2019)** forces the ordering: define the estimand (5 attributes including ICE strategy) BEFORE choosing the analysis. The missing-data strategy maps to the ICE handling strategy:

- **Treatment policy** ICE strategy + missing post-ICE data -> reference-based MI (J2R typical)
- **Hypothetical** ICE strategy -> MMRM under MAR; g-computation
- **Composite** ICE strategy -> ICE becomes part of endpoint; no missing-data problem for that subject
- **While-on-treatment** ICE strategy -> pre-ICE values only; censored at ICE
- **Principal stratum** -> latent stratum, requires Bayesian or sensitivity over unverifiable assumptions

## Missing-Data Mechanisms

| Mechanism | Definition | Testable? | Valid method |
|-----------|------------|-----------|--------------|
| MCAR | Independent of all data | Partially (Little's 1988 test) | Complete-case unbiased but loses power |
| MAR | Depends on observed data only | NOT testable | MMRM under MAR; MI under MAR |
| MNAR | Depends on unobserved values | NOT testable | Sensitivity analysis (J2R, CR, CIR, tipping point, pattern-mixture, selection model) |

**Critical philosophical point: MAR vs MNAR cannot be distinguished from observed data alone.** This is fundamental. Pre-specify the assumed mechanism in the SAP based on clinical reasoning, not the data.

## Algorithmic Taxonomy

| Method | Estimand strategy | Identification | Variance | Strength | Fails when |
|--------|-------------------|----------------|----------|----------|------------|
| Complete-case analysis | MAR or MCAR | MAR | Standard | Simple; valid under MCAR | Loses power; biased under MAR with informative covariates |
| LOCF | Implicit MNAR | Assumes flat post-ICE trajectory | Standard | Historically used | Biased even under MCAR (Mallinckrodt 2008); NRC 2010 rejects |
| MMRM with UN+KR | Hypothetical via MAR | MAR | Kenward-Roger SE | FDA-favoured continuous longitudinal | High differential dropout makes MAR implausible |
| Multiple imputation (MAR) | Hypothetical via MAR | MAR | Rubin's rules | Flexible; handles arbitrary patterns | sample_posterior must be enabled; only works with default estimator in sklearn |
| Reference-based MI (J2R, CR, CIR, LMCF) | Treatment policy / MNAR sensitivity | Clinical narrative (e.g., "after withdrawal patient resembles placebo") | Cro 2019 information-anchored vs Wolbers 2022 frequentist (active debate) | FDA-acceptable; clinician-interpretable | Variance choice contested; Rubin over-conservative, jackknife may inflate Type-I |
| Tipping-point delta-adjustment | MNAR sensitivity | Pre-specified delta function | Standard | Direct regulatory question: "how bad would missing data have to be?" | Delta interpretation depends on scale |
| Pattern-mixture with CCMV | Pattern-mixture MNAR | "Missing pattern resembles completer pattern" | Multiple imputation | Identifies missing cells via restriction | CCMV may be implausible if completers are atypical |
| Pattern-mixture with NCMV | Pattern-mixture MNAR | "Missing pattern resembles neighbouring pattern" | MI | Less extreme assumption than CCMV | Choice of "neighbouring" is ambiguous |
| Pattern-mixture with ACMV | Pattern-mixture MNAR (equivalent to MAR) | "Missing pattern resembles available cases" | MI | Equivalent to MAR (Molenberghs 1998) | Reduces to standard MAR analysis |
| Selection model (Diggle-Kenward 1994) | MNAR | Joint normal outcome + logistic dropout | Likelihood | Theoretically elegant | Conclusions driven by untestable parametric assumptions; FDA discouraged |
| Retrieved-dropout MI | Treatment policy | Sampling from observed post-discontinuation data | MI variance | Empirically grounded (no model assumption for missing) | Requires actual post-ICE data collection |

**Postdoc reading list:**

- NRC 2010 *Prevention and Treatment of Missing Data in Clinical Trials* (National Academies)
- Little RJA, D'Agostino R, Cohen ML et al 2012 *NEJM* 367:1355 (NRC summary)
- Mallinckrodt CH 2008/2014 Drug Information Journal / TIRS (MMRM case)
- Carpenter JR, Roger JH, Kenward MG 2013 *J Biopharm Stat* 23:1352 (reference-based MI)
- Cro S, Carpenter JR, Kenward MG 2019 *JRSS-A* 182:623 (information-anchored variance)
- Bartlett JW 2021 *Stat Biopharm Res* 15(1):178 (frequentist variance counter)
- Wolbers M, Noci A, Delmar P et al 2022 *Pharm Stat* 21:1262 (CMI+jackknife)
- Permutt T 2016 *Stat Med* 35:2876 (analyst as adversary; tipping point)
- Diggle PJ, Kenward MG 1994 *JRSS-C* 43:49 (selection model + canonical critique)
- Molenberghs G, Michiels B, Kenward MG, Diggle PJ 1998 *Stat Neerl* 52:153 (pattern-mixture identifying restrictions; ACMV = MAR under monotone)
- Olarte Parra C, Bartlett JW, Daniel RM 2022 *Stat Biopharm Res* (MMRM-MAR IS a hypothetical estimator)

## Decision Tree by Scenario

| Scenario | Recommended approach | Why |
|----------|---------------------|-----|
| Continuous endpoint, monotone missingness, MAR plausible | MMRM with UN + Kenward-Roger via R `mmrm` | FDA-favoured Mallinckrodt 2008/2014 default |
| Continuous endpoint, ICE = treatment discontinuation, treatment policy estimand | Hybrid: J2R for discontinuation ICEs, MMRM-MAR for other (Aprocitentan 2024 precedent) | De facto FDA standard 2024-2025 |
| Continuous endpoint, treatment policy + post-ICE data collected | Retrieved-dropout MI (Wegovy STEP precedent) | Empirically grounded; FDA 2025 obesity guidance endorses |
| Continuous endpoint, MNAR sensitivity required | J2R via `rbmi` with both Rubin's variance AND frequentist (CMI+jackknife) | Cro vs Bartlett debate; report both for safety |
| Continuous endpoint, tipping-point analysis | Permutt 2016 delta-adjustment in active arm only | Direct regulatory question; pre-specify delta range |
| Binary endpoint with missing primary | MI with logistic imputation; modified Poisson for marginal RR | Per FDA 2023 covariate adjustment |
| Time-to-event with informative censoring | IPCW (Robins) or sensitivity under composite | Censoring-as-event composite; cite Lewis 2023 |
| Very high missingness (>40%) | Report as hypothesis-generating; multiple sensitivity analyses | NRC 2010 caveat |
| Aducanumab-style differential dropout pattern | MAR primary is questionable; treatment-policy with reference-based MI primary | Lessons from FDA 6-1 AdCom decision (2021) |
| Pediatric trial with hard-to-retain population | Retrieved-dropout MI + Bayesian extrapolation from adult data | FDA 2025 obesity pediatric extension guidance |

## MMRM Under MAR -- The FDA-Preferred Continuous Analysis

**Goal:** Estimate the treatment-by-visit contrast at the primary timepoint for a continuous longitudinal endpoint under MAR with valid Type-I control in small/moderate trials.

**Approach:** Fit MMRM with unstructured covariance, REML, Kenward-Roger DF correction via the Roche/openpharma `mmrm` package; pre-specify the convergence fallback hierarchy in the SAP.

**Mallinckrodt 2008/2014:** for continuous longitudinal endpoints under monotone (or near-monotone) MAR, an MMRM with treatment + visit + treatment-by-visit + baseline + baseline-by-visit, UN covariance, REML, contrast at primary timepoint -- is consistent and FDA-preferred.

```r
library(mmrm)
fit <- mmrm(
    formula = change_from_baseline ~ baseline + arm * visit + us(visit | subject),
    data = trial_data,
    method = "Kenward-Roger-Linear",  # matches SAS PROC MIXED bit-for-bit
    reml = TRUE
)
summary(fit)
```

**Kenward-Roger flavour question:**

- `method = "Kenward-Roger"` -- full second-order KR (Kenward-Roger 1997 *Biometrics* 53:983)
- `method = "Kenward-Roger-Linear"` -- drops second-order Cholesky-derivative term; matches SAS PROC MIXED bit-for-bit
- Most regulatory submissions use Kenward-Roger-Linear for SAS-R reproducibility

**Convergence fallback hierarchy** (pre-specify in SAP):

1. UN with KR (preferred)
2. UN with Satterthwaite (if KR fails)
3. Heterogeneous Toeplitz (k+1 parameters)
4. AR(1) with heterogeneous variances
5. CS with heterogeneous variances (last resort)

**The Olarte Parra unification (2022):** MMRM under MAR IS a causal hypothetical estimand by g-formula equivalence under specific identifying assumptions. The issue is articulation, not statistical machinery.

## Reference-Based Multiple Imputation -- The rbmi Framework

**Carpenter-Roger-Kenward 2013** operationalises MNAR sensitivity as clinical narrative, not numeric delta:

- **J2R:** "after withdrawal the patient instantly resembles the placebo arm"
- **CR:** "the entire post-baseline trajectory copied from placebo, with subject's baseline deviation preserved"
- **CIR:** "the patient retains on-treatment increment but trends with placebo arm thereafter, anchored at last observed value"
- **LMCF:** "patient stays at last on-treatment mean"

```r
library(rbmi)

# Define imputation model
vars <- set_vars(
    outcome = 'CHG',  # change from baseline
    visit = 'AVISIT',
    subjid = 'USUBJID',
    group = 'ARM',
    covariates = c('BASE', 'STRATA1'),
    method = method_bayes(n_samples = 100)
)

# Draws -> Impute -> Analyse -> Pool pipeline
draws_obj <- draws(data = trial_data, vars = vars)
imputed_j2r <- impute(draws_obj,
                      references = c('Active' = 'Placebo', 'Placebo' = 'Placebo'))
analyses <- analyse(imputed_j2r,
                    fun = ancova,
                    vars = list(outcome = 'CHG', visit = 'AVISIT',
                                group = 'ARM', covariates = c('BASE')))
result <- pool(analyses)  # Rubin's rules pooling
summary(result)
```

### Inference engine choice (rbmi)

1. **Bayesian MI + Rubin's rules** (historical default) -- information-anchored variance
2. **Approximate Bayesian via REML + bootstrap** -- frequentist variance
3. **Conditional mean imputation + jackknife** (Wolbers 2022) -- deterministic; FDA-friendly
4. **BMLMI** (Lipkovich-Ratitch) -- bootstrapped MI with within/between decomposition

## The Variance Debate -- Cro vs Bartlett

**The single most active argument in current biostatistics.**

**Cro/Carpenter/Kenward 2019 *JRSS-A* 182:623:** Rubin's-rules variance applied to J2R/CR/CIR is approximately information-anchored — the relative loss of information from missingness in sensitivity analysis matches the relative loss in MAR primary analysis. True repeated-sampling variance is "information positive" because reference-based MI borrows from reference arm, reducing marginal variance of active arm BELOW what an MAR analysis with same missingness would give.

**Philosophical position (Cro et al):** a sensitivity analysis should not import information the primary analysis didn't have; if borrowing from placebo makes active CI narrower, it's no longer anchored.

**The clearer framing for postdocs:** the dispute is NOT about "which variance is conservative" — it is about "which variance answers the right question." Cro: if the sensitivity analysis is meant to assess robustness to MAR, the variance should be the one that matches the informational scope of the primary MAR analysis (Rubin's, which under-states borrowing). Bartlett: if J2R is taken as the true data-generating mechanism, then the variance of inference under that mechanism is the frequentist (jackknife) variance, which delivers nominal Type-I. Both can be correct under different framings of "what is the sensitivity analysis for?"

**Bartlett 2021 *Stat Biopharm Res* 15(1):178 + Wolbers 2022 *Pharm Stat* counter:** if J2R is the actual sampling model under which inference is made, the correct frequentist variance is the one delivering nominal Type-I and coverage — the jackknife/bootstrap variance, NOT Rubin's. Simulations: Bayesian MI + Rubin's gives Type-I 0.9-2.5% (over-conservative); CMI+jackknife gives 4.84-4.96% (nominal) under J2R.

**Regulatory practice 2024-2026:**

- EMA tolerates either
- FDA reviewers increasingly flag Rubin's variance under reference-based MI as needing a frequentist sensitivity analysis in addition
- **Safe approach:** report both — Bayesian MI + Rubin's as primary (per Cro), CMI+jackknife frequentist as supportive (per Wolbers)

## Permutt Tipping-Point Analysis

**Permutt 2016 *Stat Med* 35:2876** (Permutt was head of FDA Division of Biometrics IV at the time): the regulator's question is not "what is a reasonable MNAR adjustment?" but "how bad would the missing data have to be in the active arm to overturn the significant primary result?"

```r
library(rbmi)

# Delta-template: per-visit, per-arm, per-pattern delta
delta_grid <- seq(0, 20, by = 2)  # delta range to scan
results <- list()
for (delta in delta_grid) {
    delta_template <- delta_template(imputed, delta = delta,
                                      dlag = c(1, 1, 1, 1))  # apply to all post-ICE visits
    analyses <- analyse(imputed, delta = delta_template, fun = ancova, vars = vars)
    pooled <- pool(analyses)
    results[[as.character(delta)]] <- pooled
}

# Find minimum delta that flips p > 0.05 -> tipping delta
```

**Delta-adjustment patterns:**

- **One-arm shift (FDA preferred):** add delta to imputed values in active arm only
- **Symmetric shift:** both arms worsened by delta (probes systematic optimism)
- **Reverse shift:** placebo improved by delta (more aggressive)

**Report tipping delta in residual SD units** (FDA preference for cross-trial comparison), not raw outcome units.

## Pattern-Mixture and Selection Models

### Pattern-mixture identifying restrictions (Little 1993, Molenberghs et al 1998)

Pattern-mixture factorises joint distribution as observed-data distribution stratified by dropout pattern × pattern probability. Introduces unidentifiable parameters for unobserved cells, resolved by restrictions:

- **CCMV (Complete-Case Missing-Value):** equates missing conditional to completer pattern's conditional
- **NCMV (Neighbouring-Case Missing-Value):** uses conditional of patients with one additional measurement
- **ACMV (Available-Case Missing-Value):** weights across all patterns where relevant components observed

**Molenberghs et al 1998 proved ACMV is exactly equivalent to MAR under monotone missingness** — so pattern-mixture under ACMV is a reparameterisation of MAR analysis.

J2R/CR/CIR/LMCF are pattern-mixture models with reference-arm-based identifying restrictions.

### Selection model (Diggle-Kenward 1994)

Joints a multivariate normal response model with logistic dropout model that depends on the unobserved current value.

**Canonical critique** (Diggle-Kenward 1994 discussion; Molenberghs/Kenward/Verbeke subsequent work): selection-model MNAR conclusions are driven not by data but by parametric assumptions that are *empirically untestable* — the difference between MAR and MNAR fit is identified entirely from joint normality assumption. A non-normal response will spuriously appear MNAR.

**FDA position:** selection models are used as sensitivity, never primary. FDA Division of Biometrics has repeatedly pushed back on selection models on the grounds that "the sponsor cannot tell me in clinical English what assumption I am being asked to accept." Reference-based MI is preferred because it is clinically articulable.

## Multiple Imputation in Python -- The sklearn Caveats

```python
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
import statsmodels.formula.api as smf
import numpy as np
import pandas as pd

n_imputations = 20  # rule: m >= 100 * FMI

# CRITICAL: sample_posterior=True only works with BayesianRidge (default estimator)
imputer = IterativeImputer(max_iter=10, random_state=0, sample_posterior=True)

results = []
for i in range(n_imputations):
    imputer.set_params(random_state=i)
    imputed = pd.DataFrame(imputer.fit_transform(df[numeric_cols]), columns=numeric_cols)
    for col in ['ARM', 'sex']:
        imputed[col] = df[col].values
    model = smf.logit('outcome ~ C(ARM, Treatment(reference="Placebo")) + age', data=imputed).fit(disp=0)
    results.append({'coef': model.params.iloc[1], 'se': model.bse.iloc[1]})

# Rubin's rules
pooled_coef = np.mean([r['coef'] for r in results])
within_var = np.mean([r['se']**2 for r in results])
between_var = np.var([r['coef'] for r in results], ddof=1)
total_var = within_var + (1 + 1/n_imputations) * between_var
pooled_se = np.sqrt(total_var)
```

**Critical caveats:**

- `sample_posterior=True` SILENTLY ignored if estimator changed from BayesianRidge (e.g., to RandomForest) -> MI degenerates to single imputation
- IterativeImputer is experimental; API may change without standard deprecation
- For mixed types (binary + continuous), consider `miceforest` or R `mice`/`rbmi`
- Imputation model must include all analysis-model predictors (congeniality per Meng 1994)
- Include outcome as predictor in imputation model but exclude from imputed variables
- Never impute treatment assignment (fully determined by randomisation)

**For confirmatory regulatory work, prefer R `rbmi` or `mice` over Python sklearn** — the SAS / R precedent is stronger, the variance theory is better-developed.

## Decisive Regulatory Cases

### Aducanumab (Biogen BLA 761178, 2021)

EMERGE and ENGAGE both stopped early for futility; EMERGE high-dose positive, ENGAGE negative. MMRM-MAR primary. FDA Office of Biostatistics (Tristan Massie review) argued futility-stop-induced missingness was NOT MAR (differential ARIA-driven unblinding). 6-1 AdCom against approval; subsequent approval over-ruled OB. Textbook case showing MAR-based primary in trial with high differential missingness is regulator-divisive.

**Lesson:** when differential dropout patterns differ by arm in clinically meaningful ways, MAR is questionable; primary should be treatment-policy with reference-based MI.

### Aprocitentan (Idorsia PRECISION, FDA 2024)

Mathur 2025 *Pharm Stat* (PMC12753554) documents the negotiation. FDA pushed back on MMRM-MAR primary; accepted compromise: stratified imputation — J2R for treatment-discontinuation ICEs, MAR-MMRM for other missingness.

**Lesson:** this hybrid is now the de facto FDA standard for treatment-policy estimands. Pre-specify in SAP.

### Wegovy/Ozempic STEP trials (Wilding 2021 *NEJM*)

Retrieved-dropout MI as primary for treatment-policy estimand. Missing body weight at week 68 imputed by sampling from observed week-68 measurements among "retrieved dropouts" (patients who discontinued semaglutide but remained in follow-up). J2R-MI as supportive.

**Lesson:** RD-MI is now standard for chronic weight management. FDA 2025 obesity guidance explicitly endorses MI as primary.

## Per-Method Failure Modes

### MMRM-MAR with high differential missingness

- **Trigger:** Dropout rate or reason differs by arm (toxic active vs tolerated placebo).
- **Mechanism:** MAR requires that conditional on observed covariates, missingness is unrelated to outcome — implausible when patients drop out *because* the treatment isn't working.
- **Symptom:** Discontinuation reasons differ qualitatively; primary p sensitive to model specification.
- **Fix:** Treatment-policy estimand with retrieved-dropout MI as primary; J2R as sensitivity; tipping-point delta in active arm.

### Rubin's variance under reference-based MI

- **Trigger:** J2R/CR/CIR with only Rubin's-rules variance reported.
- **Mechanism:** Rubin's variance is information-anchored but over-conservative (Cro 2019).
- **Symptom:** Power loss vs frequentist variance.
- **Fix:** Report frequentist variance via CMI+jackknife as supportive (Wolbers 2022); cite both Cro and Bartlett.

### `sample_posterior=False` in sklearn IterativeImputer

- **Trigger:** Default behaviour ignored OR estimator changed from BayesianRidge.
- **Mechanism:** Imputations are point predictions, near-identical across draws; between-imputation variance ~ 0.
- **Symptom:** Artificially narrow Rubin's-pooled CIs.
- **Fix:** Always `sample_posterior=True`; verify estimator is BayesianRidge; consider R `rbmi`.

### Diggle-Kenward selection model as primary

- **Trigger:** Selection model used as primary MNAR analysis.
- **Mechanism:** MNAR vs MAR distinction driven by joint normality assumption, not data.
- **Symptom:** Reviewer asks "what clinical assumption am I being asked to accept?" and sponsor cannot answer in clinical English.
- **Fix:** Switch to pattern-mixture (reference-based MI) for primary; selection model as supportive sensitivity only.

### LOCF as primary or "conservative"

- **Trigger:** SAP specifies LOCF as primary analysis or as "conservative" sensitivity.
- **Mechanism:** LOCF is biased even under MCAR (Mallinckrodt 2008); discards uncertainty.
- **Symptom:** Reviewer cites NRC 2010 Rec 11 against LOCF.
- **Fix:** MMRM under MAR as primary; reference-based MI for MNAR sensitivity; LOCF only as historical comparison if at all.

### Imputation model uncongenial with analysis model

- **Trigger:** Analysis model includes treatment-by-covariate interaction but imputation model does not.
- **Mechanism:** Uncongeniality biases estimates and invalidates Rubin's variance pooling.
- **Symptom:** Pooled SE smaller than independent bootstrap suggests.
- **Fix:** Imputation model must be at least as flexible as analysis model (Meng 1994).

### Tipping delta in raw units only

- **Trigger:** Tipping-point analysis reports delta in raw outcome scale.
- **Mechanism:** Hard to compare across trials; clinical plausibility unclear.
- **Symptom:** Reviewer asks "how many SDs is that?"
- **Fix:** Report in residual SD units (FDA preference); both raw and standardised.

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| MMRM-MAR CI overlaps null; J2R MI CI excludes null | Reference-based MI borrows information from reference arm, reducing active-arm marginal variance below what MAR analysis with same missingness would give (Cro 2019) | Both valid under their assumptions; choose by clinical plausibility (MAR vs MNAR after differential dropout); report both; cite estimand strategy |
| Bayesian MI + Rubin's variance Type-I ~1%; CMI+jackknife Type-I ~5% under J2R | Rubin's information-anchored, over-conservative (Cro 2019); jackknife frequentist nominal (Wolbers 2022); active methodological debate | Report both; flag as Cro vs Bartlett debate; FDA increasingly accepts jackknife as supportive |
| LOCF "conservative" sensitivity gives smaller effect than MMRM-MAR | LOCF assumes flat post-ICE trajectory; biased even under MCAR (Mallinckrodt 2008) | Replace LOCF with reference-based MI; cite NRC 2010 Rec 11; reframe sensitivity as "MNAR robustness" not "conservative" |
| Tipping delta is small relative to MCID | MAR plausibility weak; primary result fragile to mild MNAR | Report tipping delta in residual SD AND relative to MCID; reconsider primary estimand toward treatment-policy with reference-based MI |
| Pattern-mixture with CCMV vs J2R conclusions differ | CCMV assumes missing pattern mirrors completers; J2R assumes mirror reference arm; different MNAR mechanisms | Choose based on which clinical scenario is more plausible (Carpenter-Roger 2013); report both as sensitivity range |
| Selection model (Diggle-Kenward) rejects MAR; pattern-mixture (reference-based MI) accepts MAR | Selection model MNAR conclusion driven by joint normality assumption (not data); pattern-mixture clinically articulable | Report pattern-mixture as primary sensitivity; selection model as supportive only; FDA prefers clinical articulation |
| MMRM with UN converges in arm A but not arm B | Convergence fragility in unstructured covariance with high dropout in one arm | Apply pre-specified SAP fallback (UN+KR -> UN+Satterthwaite -> heterogeneous Toeplitz -> AR(1)); document in CSR |

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| m >= 100 * FMI imputations | von Hippel 2020 *Sociol Methods Res* 49:699 | Stable pooled SE; with 40% missingness and FMI~0.3, m=30 needed |
| Missing > 40% on key variable | NRC 2010 | MI under MAR unreliable; treat as hypothesis-generating |
| Kenward-Roger DF correction for MMRM with UN | Kenward-Roger 1997 | Without it MMRM-REML under-covers; Type-I inflates 1-2 pp |
| Tipping delta in residual SD units | FDA Division of Biometrics preference | Cross-trial comparison |
| Rubin's vs frequentist for reference-based MI | Cro 2019 vs Wolbers 2022 | Report both for regulatory safety |
| Pre-specify ICE strategy in SAP | ICH E9(R1) | Estimand-before-method; cite Kahan 2023 |
| Examine DS domain for differential dropout | NRC 2010 implicit | If dropout differs by arm, MAR is suspect |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| LOCF reported as "conservative" | Persistent misconception | LOCF is biased under MCAR (Mallinckrodt 2008); cite NRC 2010 |
| MAR primary in trial with differential dropout | MAR plausibility not verified | Examine DS; if differential, switch to treatment-policy + reference-based MI |
| Selection model used as primary | Untestable parametric assumption drives MNAR conclusion | Pattern-mixture (reference-based MI) as primary; selection model as supportive |
| sample_posterior silently ignored | Estimator changed from BayesianRidge in sklearn | Verify estimator; consider R `rbmi` or `mice` |
| Imputation model missing interaction term | Uncongeniality with analysis model | Include analysis-model predictors AND interactions |
| Tipping delta in raw outcome scale only | Hard to compare | Report in residual SD units alongside raw |
| Rubin's variance for J2R without frequentist sensitivity | Cro 2019 over-conservative | CMI+jackknife frequentist as supportive (Wolbers 2022) |
| MMRM with CS forced after UN convergence failure | Pre-specification of fallback missing | Pre-specify hierarchy in SAP; document deviation if invoked |
| MI without outcome as predictor in imputation model | Underestimates association | Include outcome as predictor; exclude from imputed variables |
| MMRM in Python (statsmodels.mixedlm) treated as FDA-equivalent | Lacks Kenward-Roger | Use R `mmrm` for confirmatory; statsmodels.mixedlm only for exploratory |

## Anticipated Reviewer Pushback

| Pushback | Response |
|----------|----------|
| "Why MAR not MNAR primary?" | Examined DS domain; dropout rates and reasons symmetric across arms; cite Mallinckrodt 2008 MMRM-MAR convention |
| "Why not LOCF as sensitivity?" | LOCF biased even under MCAR; cite NRC 2010 Rec 11; J2R and tipping-point are clinically articulable alternatives |
| "Rubin's variance for J2R?" | Cite Cro 2019 information-anchored argument as rationale; report frequentist (CMI+jackknife) as supportive per Wolbers 2022 |
| "Tipping delta plausibility?" | Tipping delta = X.X in residual SD units; X.X exceeds MCID of Y.Y; deemed clinically implausible |
| "Selection model sensitivity?" | Provided as supportive; cite Diggle-Kenward 1994 critique that conclusions are driven by joint normality assumption; pattern-mixture is the primary sensitivity |
| "Imputation model congeniality?" | Imputation model includes all analysis-model predictors and interactions per Meng 1994 |
| "ICE strategy?" | Pre-specified per ICH E9(R1): treatment policy primary, hypothetical sensitivity, with explicit ICE mechanism in protocol |
| "Why m=20 imputations?" | Computed m >= 100 * FMI; observed FMI=0.20 -> m=20 sufficient (von Hippel 2020) |
| "Retrieved-dropout MI with sparse post-ICE data?" | If post-ICE retrieval is <50% complete, retrieved-dropout MI degrades to reference-based; switch to J2R as primary with reference-based MI as named sensitivity per Aprocitentan precedent |
| "Why is the reference arm chosen for J2R clinically plausible?" | Per ICH E9(R1), reference-based MI assumes post-discontinuation trajectory mirrors the placebo arm; this matches the clinical scenario of "patient stops treatment due to AE and returns to standard-of-care baseline." Documented in protocol with medical-monitor sign-off. |

## References

- Bartlett JW. 2021. Reference-based multiple imputation -- what is the right variance and how to estimate it. *Stat Biopharm Res* 15(1):178-186.
- Carpenter JR, Roger JH, Kenward MG. 2013. Analysis of longitudinal trials with protocol deviation: a framework for relevant, accessible assumptions, and inference via multiple imputation. *J Biopharm Stat* 23:1352-1371.
- Cro S, Carpenter JR, Kenward MG. 2019. Information-anchored sensitivity analysis: theory and application. *JRSS-A* 182:623-645.
- Diggle PJ, Kenward MG. 1994. Informative drop-out in longitudinal data analysis. *JRSS-C* 43:49-93.
- EMA. 2010. Guideline on Missing Data in Confirmatory Clinical Trials. EMA/CPMP/EWP/1776/99 Rev.1.
- ICH. 2019. E9(R1) Addendum on Estimands and Sensitivity Analysis.
- Kenward MG, Roger JH. 1997. Small sample inference for fixed effects from REML. *Biometrics* 53:983-997.
- Little RJA. 1988. A test of missing completely at random for multivariate data with missing values. *JASA* 83:1198-1202.
- Little RJA. 1993. Pattern-mixture models for multivariate incomplete data. *JASA* 88:125-134.
- Little RJA, D'Agostino R, Cohen ML et al. 2012. The prevention and treatment of missing data in clinical trials. *NEJM* 367:1355-1360.
- Mallinckrodt CH, Lane PW, Schnell D, Peng Y, Mancuso JP. 2008. Recommendations for the primary analysis of continuous endpoints in longitudinal clinical trials. *Drug Information Journal* 42:303-319.
- Mathur AR et al 2025. Regulatory experiences with multiple imputation. *Pharm Stat*.
- Meng XL. 1994. Multiple-imputation inferences with uncongenial sources of input. *Stat Sci* 9:538-558.
- Molenberghs G, Michiels B, Kenward MG, Diggle PJ. 1998. Monotone missing data and pattern-mixture models. *Stat Neerl* 52:153-161.
- NRC. 2010. *The Prevention and Treatment of Missing Data in Clinical Trials*. National Academies Press.
- Olarte Parra C, Bartlett JW, Daniel RM. 2022. Hypothetical estimands: a unification of causal inference and missing data methods. *Stat Biopharm Res*.
- Permutt T. 2016. Sensitivity analysis for missing data in regulatory submissions. *Stat Med* 35:2876-2879.
- von Hippel PT. 2020. How many imputations are needed (a two-stage calculation using a quadratic rule). *Sociol Methods Res* 49:699-718.
- Wolbers M, Noci A, Delmar P, Gower-Page C, Yiu S, Bartlett JW. 2022. Standard and reference-based conditional mean imputation. *Pharm Stat* 21:1262-1281.

## Related Skills

- clinical-biostatistics/trial-reporting - Estimand framework + CONSORT 2025 missing-data item 21c
- clinical-biostatistics/logistic-regression - Adjusted logistic with MI for binary endpoints
- clinical-biostatistics/effect-measures - Effect-measure CIs under MI pooling
- clinical-biostatistics/cdisc-data-handling - DS domain reasoning for missingness mechanism
- clinical-biostatistics/survival-analysis - Informative censoring as missing-data analogue
- clinical-biostatistics/adaptive-designs - Interim missing-data assumptions

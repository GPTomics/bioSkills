---
name: bio-clinical-biostatistics-logistic-regression
description: Performs logistic regression for clinical trial outcomes (binary, ordinal, multinomial) with marginal-vs-conditional estimand reporting per FDA 2023 covariate adjustment guidance, g-computation/standardisation for marginal effects, modified Poisson for RR, Brant test for proportional odds, Firth penalty for separation, and Hauck-Donner detection. Use when modeling binary or ordinal endpoints in confirmatory or exploratory clinical trials.
tool_type: python
primary_tool: statsmodels
---

## Version Compatibility

Reference examples tested with: statsmodels 0.14+, scipy 1.12+, numpy 1.26+, pandas 2.1+, firthmodels 0.3+, marginaleffects 0.0.13+ (Python) / 0.20+ (R). R packages cited: RobinCar, marginaleffects, brant, MASS, VGAM.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Logistic Regression for Clinical Outcomes

**"Model clinical outcomes with logistic regression"** -> Estimate the marginal or conditional treatment effect on a binary or ordinal endpoint using a model that respects randomisation stratification, declares its estimand, and survives covariate misspecification.

**Conditional vs marginal (the non-collapsibility subtlety):** the OR is non-collapsible. The conditional OR from logistic regression is a *different parameter* than the marginal OR, even when there is NO confounding and randomisation is perfect. This is mathematical, not statistical bias. FDA 2023 favours marginal RD (via g-computation) for primary reporting to avoid parameter ambiguity. See clinical-biostatistics/effect-measures and Permutt 2020.

## Algorithmic Taxonomy

| Approach | Estimand | Inference | Strength | Fails when |
|----------|----------|-----------|----------|------------|
| Unadjusted logistic / chi-square | Marginal OR | Wald or LR | Simple; transparent | Loses efficiency vs adjusted (Senn 2013); inflates SE under stratified randomisation (Kahan-Morris 2012) |
| Logistic with covariates (ML, Wald CI) | **Conditional** log-OR | Wald | Standard; widely available | Conditional OR != marginal OR due to non-collapsibility (Permutt 2020); not the FDA 2023 primary estimand |
| Logistic + g-computation / standardisation | **Marginal** RD/RR/OR | Influence function or bootstrap SE | FDA 2023 recommended primary estimand for binary | Requires correct outcome model AND post-fit standardisation; needs robust SE machinery |
| Targeted Maximum Likelihood (TMLE) | Marginal RD/RR/OR | Influence function | Provably efficient; doubly robust in observational | Implementation heavier; mostly R (`tmle`, `tmle3`); rare in confirmatory submissions |
| Modified Poisson with sandwich SE | Marginal RR | HC1/HC3 sandwich | Direct RR estimation when prevalence >10% | Slightly less efficient than log-binomial when log-binomial converges |
| Log-binomial regression | Marginal RR | Wald | Direct RR estimation | Frequent convergence failure when predicted risk near 1 |
| Firth penalised logistic | Conditional OR (penalised) | Penalised LR test preferred | Handles separation, rare events (<5% prevalence) | Wald CI/p liberal; must use PLR test (Heinze-Schemper 2002) |
| Ordinal logistic (proportional odds) | Common conditional OR across cut-points | Wald or LR | Preserves ordering information | Proportional odds assumption violation (Brant test) |
| Partial proportional odds | PO holds for some covariates, not others | Hybrid | Salvages ordinal model when PO fails for one predictor | Increased complexity; harder interpretation |
| Multinomial logistic | Per-category log-OR | Wald | No PO assumption needed | Loses efficiency; harder communication |

**Postdoc reading list:** Permutt 2020 *Stat Biopharm Res* 12:45 (conditional vs marginal estimand); FDA May 2023 Final Guidance "Adjusting for Covariates in RCTs" (the regulatory rulebook); Tsiatis et al 2008 *Stat Med* 27:4658 (robust ANCOVA framework); Senn 2013 *Stat Med* 32:1439 (precision from adjustment even under perfect balance); Kahan-Morris 2012 *Stat Med* 31:328 (must adjust for stratification factors).

## Decision Tree by Scenario

| Scenario | Recommended approach | Why |
|----------|---------------------|-----|
| RCT, binary primary endpoint, no covariate adjustment in SAP | Unadjusted logistic with marginal OR; chi-square supportive | Simple regulatory case; consider whether covariate adjustment would gain power (Senn 2013) |
| RCT, binary primary, covariates pre-specified, FDA 2023 compliant | **Logistic adjusted + g-computation for marginal RD**; conditional OR supportive | FDA 2023 final: marginal estimand for primary; mention HC3 sandwich SE |
| RCT, stratified randomisation (sex/region/severity) | Include strata as covariates in logistic; the analytic decision is non-optional | Kahan-Morris 2012: ignoring strata inflates Type-I error up to 30% |
| Common outcome (prevalence >10%) where RR is the policy quantity | Modified Poisson + HC1/HC3 OR log-binomial regression | Avoid OR misinterpretation; Zou 2004; cite EMA 2015 on covariate adjustment |
| Rare event (<5% prevalence) or separation observed | Firth penalty + penalised LR test for p-value | Heinze-Schemper 2002; Wald p from Firth is liberal |
| Ordinal outcome, PO assumption supportable | Proportional odds; verify with Brant test (R `brant::brant`) | Most efficient when PO holds; common in toxicity grading and PROs |
| Ordinal outcome, PO fails on one or two predictors | Partial PO model (R `VGAM::vglm(..., cumulative(parallel=FALSE~X))`) | Salvages most efficiency; only the offending coefficient gets per-cut-point estimate |
| Ordinal outcome, PO fails widely | Multinomial logistic | Worst case; loses ordering info but valid |
| Single arm observational with strong confounding | Logistic adjusted + g-computation + propensity weighting | Doubly robust; cite Hernan-Robins; consider TMLE |
| Binary endpoint, longitudinal repeated measures | GEE or generalised linear mixed model | Sandwich SE for GEE; mixed-model OR is subject-specific not marginal |

## Standard Workflow

**Goal:** Fit a covariate-adjusted logistic regression for a binary clinical endpoint with explicit reference category and pre-specified covariates.

**Approach:** Use the formula API for automatic categorical handling and explicit reference; report both the conditional log-OR and the marginal RD via g-computation.

```python
import statsmodels.formula.api as smf
import pandas as pd
import numpy as np

# CRITICAL: set reference category explicitly. statsmodels default is alphabetical,
# so 'Active' sorts before 'Placebo' and the OR direction silently flips.
model = smf.logit(
    'outcome ~ C(ARM, Treatment(reference="Placebo")) + age + C(sex) + baseline_score',
    data=df
).fit()

# Conditional ORs (Wald)
or_table = pd.DataFrame({
    'OR': np.exp(model.params),
    'Lower_CI': np.exp(model.conf_int()[0]),
    'Upper_CI': np.exp(model.conf_int()[1]),
    'p_value': model.pvalues
})

# Marginal RD via g-computation (the FDA 2023-recommended primary estimand):
df_active = df.assign(ARM='Active')
df_placebo = df.assign(ARM='Placebo')
risk_active = model.predict(df_active).mean()
risk_placebo = model.predict(df_placebo).mean()
marginal_rd = risk_active - risk_placebo
# For SE, use influence-function or bootstrap; see RobinCar / marginaleffects in R
```

**The reference-category trap is the single most common silent bug.** `smf.logit('y ~ C(ARM)')` uses alphabetical ordering, so `ARM='Active'` becomes 1 and `ARM='Placebo'` becomes 0 -- but with `ARM=['Active','Placebo','Control']`, 'Active' is the reference and 'Placebo' is the comparison. Always pass `Treatment(reference="Placebo")` or equivalent. R `glm(y ~ relevel(ARM, ref='Placebo'))`.

## Marginal vs Conditional Estimand -- The FDA 2023 Pivot

**The single most important methodological shift in clinical biostatistics 2020-2025.** Permutt 2020 established and FDA 2023 codified: the maximum-likelihood coefficient on Z in `glm(Y ~ Z + X, family=binomial)` is the **conditional log-OR** -- the OR comparing Z=1 to Z=0 holding X fixed. Due to OR non-collapsibility, this is a *different parameter* than the **marginal log-OR** (the OR comparing the entire treated population to the entire control population), even under perfect randomisation.

**The FDA 2023 final guidance** ("Adjusting for Covariates in Randomized Clinical Trials," May 2023) endorses covariate-adjusted nonlinear-model analysis *provided the analyst targets a marginal estimand*. Reporting the conditional OR from a multivariable logistic regression as "the treatment effect" silently changes the estimand from what was pre-specified.

### G-computation / Standardisation

```python
# Marginal effects on three scales via g-computation:
df_z1 = df.assign(ARM='Active')
df_z0 = df.assign(ARM='Placebo')
p_z1 = model.predict(df_z1)
p_z0 = model.predict(df_z0)

marg_p1 = p_z1.mean()
marg_p0 = p_z0.mean()
marg_rd = marg_p1 - marg_p0
marg_rr = marg_p1 / marg_p0
marg_or = (marg_p1 / (1 - marg_p1)) / (marg_p0 / (1 - marg_p0))
```

**For valid SE:** the delta-method/influence-function variance and HC3 sandwich are required. Python `marginaleffects` v0.0.13+ implements this via `marginaleffects.avg_comparisons(model, variables='ARM', vcov='HC3')`. R is more mature: `marginaleffects::avg_comparisons` (Arel-Bundock-Greifer-Heiss 2024 *JSS* 111:9), `RobinCar` (purpose-built for FDA 2023), `riskCommunicator`.

**Tsiatis et al 2008 robustness guarantee:** under randomisation Z ⊥ X, the g-computation marginal estimator is *consistent* for the marginal ATE even if the outcome model is misspecified. Efficiency depends on model quality; consistency does not. This is why FDA 2023 accepts the marginal RD via g-computation without requiring proof of correct logistic mean structure.

**Reporting template (post-FDA-2023):**

> "Primary estimand: marginal risk difference of -8.5 percentage points (95% CI -12.7 to -4.2, HC3 SE), computed by g-computation/standardisation from a logistic regression adjusted for age, sex, and baseline severity. Supportive: conditional OR 0.48 (95% CI 0.34-0.67, Wald) from the same model."

## Modified Poisson for Common Outcomes -- Direct RR

When prevalence > 10% and the policy quantity is RR (not OR), modified Poisson with sandwich SE is the de facto modern standard (Zou 2004 *AJE* 159:702):

```python
import statsmodels.api as sm
import numpy as np

X = sm.add_constant(df[['treatment', 'age', 'baseline_score']])
poisson_model = sm.GLM(df['outcome'], X, family=sm.families.Poisson()).fit(cov_type='HC1')
rr = np.exp(poisson_model.params)
rr_ci = np.exp(poisson_model.conf_int())
```

`cov_type='HC1'` corrects the over-dispersion that Poisson inherently assumes. Without it the variance is wrong because binary data are NOT Poisson. For n < 250, HC3 is preferred (Long-Ervin 2000). **What postdocs argue about:** modified Poisson vs log-binomial -- log-binomial directly models RR but frequently fails to converge when predicted risk approaches 1; modified Poisson always converges but is marginally less efficient when log-binomial works.

## Proportional Odds and Brant Test

**Proportional odds (PO)** assumption: the effect of each predictor is constant across all cut-points of the ordinal outcome. Must be tested or the model parameters are invalid.

```python
from statsmodels.miscmodels.ordinal_model import OrderedModel
from statsmodels.api import MNLogit
import pandas as pd

df['severity'] = pd.Categorical(df['severity'], categories=['mild', 'moderate', 'severe'], ordered=True)

po_model = OrderedModel.from_formula('severity ~ treatment + age', data=df, distr='logit').fit(method='bfgs', disp=0)
mn_model = MNLogit.from_formula('severity ~ treatment + age', data=df).fit(disp=0)

# LR test for PO vs MN (omnibus)
lr_stat = 2 * (mn_model.llf - po_model.llf)
lr_df = mn_model.df_model - po_model.df_model
from scipy.stats import chi2
lr_p = 1 - chi2.cdf(lr_stat, lr_df)

# Per-coefficient Brant test in R: brant::brant(po_model_object)
# The Brant test localises which predictor(s) violate PO -- essential before deciding on partial PO
```

**The omnibus LR test is necessary but not sufficient.** Brant test (Brant 1990 *Biometrics* 46:1171; R `brant::brant`) provides per-coefficient PO tests, identifying *which* predictor violates PO. If only one or two predictors violate PO, fit a **partial proportional odds model** (R `VGAM::vglm(..., cumulative(parallel=FALSE~X1+X2))`) rather than abandoning the model entirely.

**OrderedModel intercept gotcha:** do NOT add an intercept term. Threshold parameters (cut-points between ordinal levels) replace the intercept. An explicit constant causes non-identifiability and optimiser failure.

## Separation and Firth Penalty

**Detection:**

```python
from firthmodels import FirthLogisticRegression, detect_separation
import numpy as np

sep_result = detect_separation(X, y)
if sep_result.separation:
    print(sep_result.summary())
# Manual signs: coefficient > 10, SE > 100, convergence warnings
```

**Firth penalty (Firth 1993 *Biometrika* 80:27):**

```python
firth = FirthLogisticRegression()
firth.fit(X, y)
or_firth = np.exp(firth.coef_)

# IMPORTANT: firth.pvalues_ is Wald, which is LIBERAL with Firth.
# Prefer penalised LRT:
# - firth.pvalues_lrt_ if available (firthmodels >=0.3)
# - Otherwise compute manually by comparing penalised log-likelihoods of nested models
```

**Heinze-Schemper 2002 *Stat Med* 21:2409** showed Wald inference from Firth is **liberal** (anti-conservative); the penalised likelihood ratio test (PLRT) is the recommended inference. The `firthmodels` package added `pvalues_lrt_` in version 0.3; in older versions, compute manually:

```python
def penalised_lrt(firth_full, firth_reduced):
    # Compute 2 * (penalised log-lik full - penalised log-lik reduced) ~ chi-square_df
    pass  # implementation depends on package version; see Heinze-Schemper 2002
```

Firth's method was originally designed for finite-sample bias reduction, not separation per se -- it adds the Jeffreys prior penalty to the likelihood, keeping coefficients finite under separation as a side effect. Also recommended for rare events (<5% prevalence) where ML bias is non-negligible.

## Hauck-Donner Effect Detection

**The Hauck-Donner effect (1977 *JASA* 72:851; revived by Yee 2022 *JASA* 117:1763):** the Wald test statistic is non-monotonic in the parameter estimate near the boundary. A large log-OR can produce a small Wald chi-square -- so the Wald test fails to reject when LR/profile-likelihood would. Common in small samples with strong predictors.

```python
# In R, detect with VGAM::hdeff() and replace Wald with profile likelihood:
# MASS::confint.glm(model) returns profile-likelihood CIs as default
# Python equivalent: bootstrap or manual profile likelihood
```

**When to suspect Hauck-Donner:** large coefficient magnitude with non-significant Wald p; large SE with finite estimate; switching from Wald to LR test changes significance. The fix is always: switch to profile likelihood or LR inference.

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| Conditional OR from logistic vs marginal RD from g-computation give different conclusions | Non-collapsibility (OR is non-collapsible; RD is collapsible) | Report marginal RD as primary per FDA 2023; conditional OR as supportive with explicit parameter label; cite Permutt 2020 |
| ML logistic diverges (large coefficient, huge SE); Firth penalised converges | Complete or quasi-complete separation in covariate-outcome relationship | Firth penalty with penalised LR test (NOT Wald p-values); cite Heinze-Schemper 2002 |
| Modified Poisson and log-binomial give different RR estimates | Log-binomial convergence fragility when predicted risk approaches 1 | Modified Poisson with sandwich SE preferred for robust convergence; cite Zou 2004 |
| Proportional-odds (cumulative logit) and multinomial give different inferences | PO assumption violated for one or more predictors (Brant test rejects) | Localise via Brant test; if 1-2 predictors violate, partial-PO model in VGAM; if widely, multinomial |
| Adjusted vs unadjusted treatment effect differ substantially | Confounding (in observational) OR non-collapsibility (in RCT) | In RCT, non-collapsibility expected (Permutt 2020); in observational, investigate confounding via DAG |
| Large OR with non-significant Wald p-value | Hauck-Donner effect: Wald non-monotone near boundary (Yee 2022) | Use profile-likelihood CI (R `MASS::confint.glm`); detect via `VGAM::hdeff()` |
| Stratified analysis significant; unstratified not | Achieved SE smaller in stratified analysis | Include stratification factors in primary model (Kahan-Morris 2012); ignoring inflates Type-I |
| Significant treatment effect on conditional model vanishes when mediator included | Adjusting for post-treatment variable on causal pathway | Never adjust for mediators in primary; use causal DAG to distinguish confounder vs mediator |

## Per-Method Failure Modes

### Reference-category silent reversal

- **Trigger:** Treatment variable coded as character with `Active` and `Placebo`; reference not explicitly set.
- **Mechanism:** statsmodels defaults to alphabetical ordering; 'Active' becomes the reference.
- **Symptom:** OR for "ARM[T.Placebo]" appears in output; direction is the opposite of intended.
- **Fix:** Always pass `C(ARM, Treatment(reference="Placebo"))` or numeric encoding with explicit comment.

### Adjusting for a mediator

- **Trigger:** Including a post-randomisation variable that lies on the causal pathway.
- **Mechanism:** Adjustment blocks the causal effect, attenuating treatment effect estimate toward null.
- **Symptom:** Significant unadjusted effect becomes non-significant after adjustment for "mechanism" variable (e.g., inflammation marker).
- **Fix:** Use causal DAG to distinguish confounder from mediator; never adjust for post-treatment variables in the primary analysis. See ICH E9(R1).

### Hauck-Donner non-monotonicity

- **Trigger:** Small samples, strong predictor, cell counts approaching zero.
- **Mechanism:** Wald test statistic is non-monotone in the parameter near the boundary.
- **Symptom:** Large OR magnitude with non-significant Wald p; switching to LR test changes significance.
- **Fix:** Use profile likelihood (R `MASS::confint.glm`) or LR inference; detect with `VGAM::hdeff()`.

### Complete or quasi-complete separation

- **Trigger:** A predictor perfectly predicts outcome in a subset of the data.
- **Mechanism:** Likelihood is monotone in that coefficient; ML estimate diverges to infinity.
- **Symptom:** Coefficient >10, SE >100, convergence warning; convergence message says "fitted probabilities numerically 0 or 1."
- **Fix:** Firth penalty (`firthmodels`); inference via penalised LR test, not Wald.

### Proportional odds violation undetected

- **Trigger:** Ordinal outcome fit with cumulative-logit (PO) model without testing assumption.
- **Mechanism:** PO assumes a constant effect across cut-points; violation invalidates model parameters.
- **Symptom:** Brant test rejects PO for one or more predictors; per-cut-point effects in a saturated model diverge.
- **Fix:** Brant test first; if PO fails on one predictor, partial PO model; if widely fails, multinomial logistic.

### Stratified randomisation not in analysis

- **Trigger:** Stratification by site/region/severity at randomisation; analysis ignores strata.
- **Mechanism:** Achieved SE smaller than calculated SE because randomisation removed between-stratum variability.
- **Symptom:** Anti-conservative Type-I error; Kahan-Morris 2012 documents up to 30% inflation in published trials.
- **Fix:** Include stratification factors as covariates in the logistic; this is non-optional per ICH E9, FDA 2023, EMA 2015.

## Model Diagnostics

| Diagnostic | Method | Threshold | Caveat |
|-----------|--------|-----------|--------|
| Discrimination | ROC-AUC (`sklearn.metrics.roc_auc_score`) | >0.7 acceptable, >0.8 good | Trial-level, not patient-individual; depends on outcome prevalence |
| Calibration -- primary | Calibration plot (observed vs predicted in deciles) | Curve along the diagonal | Visual; preferred over H-L |
| Calibration -- secondary | Hosmer-Lemeshow chi-square | p > 0.05 | Low power n<200; oversensitive n>2000 |
| Pseudo R-squared | model.prsquared (McFadden) | >0.2 excellent | NOT comparable to OLS R-squared |
| Events per variable | n_events / n_covariates | >=10 EPV (Peduzzi 1996) | Below this: bias, overfitting; consider Firth |

### Hosmer-Lemeshow

```python
from scipy.stats import chi2
import pandas as pd

def hosmer_lemeshow(y_true, y_pred, n_groups=10):
    df_hl = pd.DataFrame({'y': y_true, 'prob': y_pred})
    df_hl['group'] = pd.qcut(df_hl['prob'], n_groups, duplicates='drop')
    grouped = df_hl.groupby('group').agg(obs=('y', 'sum'), n=('y', 'count'), pred=('prob', 'mean'))
    grouped['expected'] = grouped['n'] * grouped['pred']
    hl_stat = (((grouped['obs'] - grouped['expected']) ** 2) /
               (grouped['n'] * grouped['pred'] * (1 - grouped['pred']))).sum()
    actual_groups = len(grouped)
    return hl_stat, 1 - chi2.cdf(hl_stat, actual_groups - 2)
```

**H-L is supplementary, not primary:** Hosmer-Lemeshow has low power n<200 (rarely rejects even for poor calibration) and oversensitive n>2000 (rejects for trivial miscalibration). The decile choice is also arbitrary -- `pd.qcut(..., duplicates='drop')` can reduce the actual number of groups when probabilities tie, changing df. Use calibration plots (smoothed observed vs predicted) as primary; Hosmer-Lemeshow as supplementary p-value.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| 10 events per variable (EPV) | Peduzzi et al 1996 *J Clin Epidemiol* 49:1373 | Below this, coefficient bias and CI mis-coverage |
| Prevalence > 10% -> prefer marginal RR via modified Poisson | Zou 2004 *AJE* 159:702 | OR overstates RR; modified Poisson directly estimates RR |
| Marginal estimand for primary regulatory analysis | FDA 2023 Final Guidance | Conditional OR is a different parameter than marginal (Permutt 2020) |
| Brant test for PO assumption | Brant 1990 *Biometrics* 46:1171 | Omnibus LR test misses per-coefficient violations |
| Firth penalty with penalised LR test, not Wald | Heinze-Schemper 2002 *Stat Med* 21:2409 | Wald is liberal under Firth penalty |
| HC3 sandwich SE for n <=250 | Long-Ervin 2000 *Am Stat* 54:217 | HC1 (Stata default) anti-conservative in small samples |
| Stratification factors in analysis when used in randomisation | Kahan-Morris 2012 *Stat Med* 31:328 | Ignoring inflates Type-I error up to 30% |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| OR direction opposite of expected | Alphabetical reference; 'Active' sorts before 'Placebo' | `C(ARM, Treatment(reference="Placebo"))` always |
| Conditional and marginal OR differ substantially | Non-collapsibility, not confounding | Cite Permutt 2020; report both with explicit labels |
| Treatment effect vanishes after adjustment | Adjusted for a mediator (post-treatment variable) | Causal DAG check; never adjust for post-treatment in primary |
| Huge coefficient, huge SE, non-significant Wald | Separation OR Hauck-Donner | Detect separation: `firthmodels.detect_separation`; switch to Firth + PLRT |
| H-L p > 0.5 with obvious miscalibration | Low power at n<200 | Use calibration plot as primary; H-L supplementary |
| H-L p < 0.001 with great-looking plot | Oversensitive at n>2000 | Use calibration plot; cite Steyerberg 2019 for cautions |
| `OrderedModel` fails to converge with intercept | Threshold parameters replace intercept | Drop the explicit constant |
| `firth.pvalues_` looks too low | Wald p from Firth is liberal | Use `firth.pvalues_lrt_` (>=0.3) or compute PLRT manually |
| GLM(family=Binomial) and Logit give same point but different output | `sm.Logit` provides `pred_table()`, `prsquared`; `sm.GLM` does not | Use Logit unless changing link functions |
| Robust SE not reported for modified Poisson | Missing `cov_type='HC1'` or 'HC3' | Always specify sandwich SE for modified Poisson |

## Anticipated Reviewer Pushback

| Pushback | Response |
|----------|----------|
| "Why is the marginal RD different from the conditional OR's implied RD?" | Non-collapsibility of OR; cite Permutt 2020 and FDA 2023. Marginal RD via g-computation is the primary per FDA 2023. |
| "Adjustment for stratification factors?" | Pre-specified strata included in the model. Cite Kahan-Morris 2012 and ICH E9. |
| "How were covariates chosen?" | Pre-specified in the SAP based on prior literature/clinical knowledge. No data-driven selection (which would inflate Type-I). |
| "Why not use stepwise selection?" | Stepwise leaks signal from the outcome and inflates Type-I. Pre-specification is the regulatory norm; cite FDA 2023. |
| "PH check for the longitudinal binary?" | Binary doesn't have PH; if longitudinal, used GEE with sandwich SE or GLMM; see ICH E9 estimand for treatment policy. |
| "Why Firth not standard logistic?" | Separation detected (or rare events <5%); standard ML diverges; Firth penalty + PLRT recommended (Heinze-Schemper 2002). |
| "Brant test result?" | PO holds for predictors X1, X2, X3; PO fails for X4 -> partial PO model fit per VGAM. |
| "Calibration?" | Calibration plot in supplement; H-L p reported as supplementary, not primary. |

## References

- Arel-Bundock V, Greifer N, Heiss A. 2024. How to interpret statistical models using marginaleffects in R and Python. *J Stat Softw* 111:9.
- Brant R. 1990. Assessing proportionality in the proportional odds model for ordinal logistic regression. *Biometrics* 46:1171-1178.
- FDA. 2023. Adjusting for Covariates in Randomized Clinical Trials for Drugs and Biological Products. Final Guidance, May 2023.
- Firth D. 1993. Bias reduction of maximum likelihood estimates. *Biometrika* 80:27-38.
- Hauck WW, Donner A. 1977. Wald's test as applied to hypotheses in logit analysis. *JASA* 72:851-853.
- Heinze G, Schemper M. 2002. A solution to the problem of separation in logistic regression. *Stat Med* 21:2409-2419.
- Kahan BC, Morris TP. 2012. Improper analysis of trials randomised using stratified blocks or minimisation. *Stat Med* 31:328-340.
- Lin DY, Wei LJ. 1989. The robust inference for the Cox proportional hazards model. *JASA* 84:1074-1078.
- Long JS, Ervin LH. 2000. Using heteroscedasticity consistent standard errors in the linear regression model. *Am Stat* 54:217-224.
- Moore KL, van der Laan MJ. 2009. Covariate adjustment in randomized trials with binary outcomes: TMLE. *Stat Med* 28:39-64.
- Peduzzi P, Concato J, Kemper E, Holford TR, Feinstein AR. 1996. A simulation study of the number of events per variable in logistic regression analysis. *J Clin Epidemiol* 49:1373-1379.
- Permutt T. 2020. Do covariates change the estimand? *Stat Biopharm Res* 12:45-53.
- Senn S. 2013. Seven myths of randomisation in clinical trials. *Stat Med* 32:1439-1450.
- Steyerberg EW. 2019. *Clinical Prediction Models* (2nd ed). Springer.
- Tsiatis AA, Davidian M, Zhang M, Lu X. 2008. Covariate adjustment for two-sample treatment comparisons in randomized clinical trials. *Stat Med* 27:4658-4677.
- Wang B, Susukida R, Mojtabai R, Amin-Esmaeili M, Rosenblum M. 2021. Model-robust inference for clinical trials that improve precision by stratified randomization and covariate adjustment. *JASA* 116:1856-1870.
- Yee TW. 2022. On the Hauck-Donner effect in Wald tests. *JASA* 117:1763-1774.
- Zou G. 2004. A modified Poisson regression approach to prospective studies with binary data. *AJE* 159:702-706.

## Related Skills

- clinical-biostatistics/cdisc-data-handling - Prepare analysis datasets from CDISC SDTM/ADaM domains
- clinical-biostatistics/effect-measures - Modern CI methods; marginal vs conditional in depth
- clinical-biostatistics/categorical-tests - Chi-square and Fisher alternatives for unadjusted tests
- clinical-biostatistics/subgroup-analysis - Interaction terms and HTE detection
- clinical-biostatistics/survival-analysis - Cox regression for time-to-event analogues
- clinical-biostatistics/multiplicity-graphical - Bretz-Maurer graphs for multiple endpoints from one logistic model
- clinical-biostatistics/trial-reporting - CONSORT 2025 and ICH E9(R1) reporting of logistic analyses

---
name: bio-clinical-biostatistics-effect-measures
description: Computes and interprets treatment effect measures (OR, RR, RD, HR, NNT) with calibrated confidence intervals (Wilson, Newcombe, Miettinen-Nurminen, MOVER, profile likelihood, Bender NNT) and reports marginal vs conditional estimands per FDA 2023 covariate adjustment guidance. Use when reporting treatment effects in confirmatory trials, comparing effect sizes across studies, or constructing forest plots.
tool_type: python
primary_tool: statsmodels
---

## Version Compatibility

Reference examples tested with: statsmodels 0.14+, numpy 1.26+, pandas 2.1+, matplotlib 3.8+, marginaleffects (Python) 0.0.13+ / (R) 0.20+. R packages cited: ratesci, exact2x2, marginaleffects, riskCommunicator, RobinCar.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Treatment Effect Measures for Clinical Trials

**"Compute treatment effect sizes"** -> Estimate the population-level treatment contrast (OR, RR, RD, HR, NNT) with a confidence interval calibrated to sample size and a clear declaration of whether the estimand is marginal or conditional under ICH E9(R1).

## Algorithmic Taxonomy

| Measure | Scale | Collapsible? | Best CI method | When to use | Fails when |
|---------|-------|--------------|----------------|-------------|------------|
| OR | Log-odds ratio | NO (non-collapsible) | Profile likelihood; Wald acceptable for n>100 per arm | Case-control (only valid measure); logistic regression default | Outcome prevalence > 10% (OR overstates RR); Hauck-Donner pathology near boundary |
| RR | Log-risk ratio | YES | Miettinen-Nurminen score; MOVER-R | Cohort, RCT with common outcomes | Sparse strata; one or both p near 0 (Wald log-RR breaks) |
| RD (absolute risk difference) | Linear probability | YES | Newcombe-Wilson hybrid; Miettinen-Nurminen | Clinically interpretable absolute scale; FDA-preferred for binary | Predictions outside [0,1] from linear models |
| HR (hazard ratio) | Log-hazard ratio | NO | Wald with profile likelihood for small n | Time-to-event with PH | PH violation (see clinical-biostatistics/survival-analysis) |
| NNT/NNH | 1/RD | n/a (derived from RD) | Bender 2001 *CCT* 22:102 (Altman 1998 base) | Communicating absolute benefit to clinicians | RD CI crosses zero (NNT becomes NNTB-infinity-NNTH) |
| Difference in RMST | Time scale | YES | Wald with delta method; pseudo-observation regression | Time-to-event with PH violation | Different max follow-up across arms (truncation tau ambiguous) |

**Postdoc reading:** Permutt 2020 *Stat Biopharm Res* 12:45 ("Do covariates change the estimand?") established that the conditional OR from a multivariable logistic regression is a *different parameter* than the marginal OR -- and the FDA May 2023 Final Guidance "Adjusting for Covariates in RCTs" requires the **marginal** estimand for primary reporting (see clinical-biostatistics/logistic-regression for g-computation/standardisation machinery).

## Decision Tree by Scenario

| Scenario | Recommended estimand + CI | Why |
|----------|---------------------------|-----|
| RCT, binary outcome, prevalence <10% | OR with Wald CI; report as primary | OR ~= RR at low prevalence; standard regulatory currency |
| RCT, binary outcome, prevalence >=10% | RR or RD via modified Poisson with HC1/HC3 sandwich SE; OR as secondary | OR substantially overstates RR; modified Poisson directly estimates RR (Zou 2004 *AJE* 159:702) |
| RCT, binary outcome, primary endpoint per FDA 2023 | Marginal RD via g-computation with sandwich SE; conditional OR as supportive | FDA 2023 final guidance: marginal estimand for primary; cite Permutt 2020 |
| Case-control study | OR only (RR unidentifiable); profile likelihood CI | OR is the only measure estimable from case-control design |
| Noninferiority on absolute scale | Miettinen-Nurminen score CI for RD | Regulatory standard for RD CIs; consistent with Pearson chi-square |
| Noninferiority on relative scale | Miettinen-Nurminen score CI for RR; Koopman 1984 acceptable | Wald log-RR has poor coverage near boundary; MN is the regulatory expectation |
| Stratified design with site/region strata | MH pooled OR with stratified MN CI (R `ratesci::scoreci(..., stratified=TRUE)`) | Preserves stratification; cite Kahan-Morris 2012 |
| Reporting NNT for clinicians | Bender 2001 *CCT* 22:102 method; report as NNTB(lower)..infinity..NNTH(upper) when CI crosses zero | Standard in BMJ/Lancet/Cochrane |
| Time-to-event with PH violation | RMST difference; cite Royston-Parmar 2013, Uno 2014 | HR is a misleading single-number summary under non-PH (see survival-analysis) |

## Crude Effect Measures from 2x2 Tables

**Goal:** Compute unadjusted OR, RR, RD from a contingency table with calibrated CIs.

**Approach:** Use Table2x2 for OR/RR with Wald CIs, but switch to score-based methods (Miettinen-Nurminen, Newcombe-Wilson) for regulatory contexts.

```python
from statsmodels.stats.contingency_tables import Table2x2
from statsmodels.stats.proportion import confint_proportions_2indep
import numpy as np

# Table layout: [[treated_event, treated_no_event], [control_event, control_no_event]]
table = np.array([[a, b], [c, d]])
t = Table2x2(table)
print('OR:', t.oddsratio, t.oddsratio_confint())        # Wald log-OR CI
print('RR:', t.riskratio, t.riskratio_confint())         # Wald log-RR CI

# Modern CI for RD (regulatory preferred):
ci_rd = confint_proportions_2indep(a, a+b, c, c+d, method='newcomb', alpha=0.05)
# 'newcomb' = Newcombe-Wilson hybrid; 'agresti-caffo' is the +1+1 +1+1 adjustment
# For Miettinen-Nurminen stratified CIs, use R `ratesci::scoreci`
```

**Table orientation is critical:** Table2x2 interprets the first column as "outcome present." `pd.crosstab` orders columns alphabetically; if outcome is coded 0/1, the table will have 0 first and the OR will be the *reciprocal* of intended. Always reorder: `cross = cross[[1, 0]]` or `cross = cross[['Yes', 'No']]`. This is a silent direction-reversing error.

**For the underlying significance test of a 2x2:** Pearson chi-square (no Yates) at n>=40 with adequate expected counts; otherwise Boschloo's exact (uniformly more powerful than Fisher's exact at the same Type-I per Mehta-Senchaudhuri 2003; Lydersen-Fagerland-Laake 2009). See clinical-biostatistics/categorical-tests for the algorithmic taxonomy.

## Modern Confidence Intervals -- the Postdoc Toolkit

Wald CIs are the textbook default but have well-documented coverage failures (Brown-Cai-DasGupta 2001 *Stat Sci* 16:101 -- "chaotic" coverage near 0 and 1 even with moderate n; Newcombe 1998a *Stat Med* 17:873 documents 11 alternatives for the difference of two proportions).

### Single proportion

| Method | When | Status |
|--------|------|--------|
| Wald | Never preferred | Defunct for serious work; coverage can drop to 0.0 (Brown-Cai-DasGupta) |
| Wilson (score) | General default | Brown-Cai-DasGupta recommended; matches Pearson chi-square test |
| Jeffreys (Beta(1/2,1/2)) | Small n | Equal-tailed Jeffreys; Bayesian with reference prior; recommended for small n |
| Clopper-Pearson | When exact guarantee required | Over-covers by 1-4 percentage points; sometimes regulatory-required |
| Agresti-Coull | Simple Wald replacement | +2/+2 adjustment; teaching default |
| Mid-p Clopper-Pearson | When CP over-coverage hurts NI margin | Less conservative than CP; slight under-coverage |

```python
from statsmodels.stats.proportion import proportion_confint
ci = proportion_confint(45, 60, alpha=0.05, method='wilson')
# Other methods: 'jeffreys', 'agresti_coull', 'beta' (Clopper-Pearson), 'normal' (Wald)
```

### Difference of two proportions (RD)

| Method | Citation | Verdict |
|--------|----------|---------|
| Wald | textbook | Poor; can produce limits outside [-1,1] |
| Newcombe-Wilson hybrid / MOVER | Newcombe 1998a | Recommended; balanced of computational simplicity and coverage |
| Agresti-Caffo (+1+1, +1+1) | Agresti-Caffo 2000 *Am Stat* 54:280 | Add 1 success + 1 failure per arm, then Wald; surprisingly good for small n |
| Miettinen-Nurminen score | Miettinen-Nurminen 1985 *Stat Med* 4:213 | **Regulatory standard for RD**; consistent with Pearson chi-square; SAS PROC FREQ `riskdiff(method=mn)` |
| Chan-Zhang exact unconditional | Chan-Zhang 1999 *Biometrics* 55:1202 | Guaranteed coverage; needed when MN under-covers near extreme p |

### Ratio of two proportions (RR)

| Method | Citation | Verdict |
|--------|----------|---------|
| Wald on log(RR) ("Katz log") | textbook (Katz 1978) | Defunct for serious work; biased when either p small |
| Koopman 1984 score | *Biometrics* 40:513 | Iterative; original score CI for RR |
| Miettinen-Nurminen score | Miettinen-Nurminen 1985 | Regulatory standard; R `ratesci::scoreci(contrast='RR')` |
| MOVER-R | Donner-Zou 2012 *Stat Methods Med Res* 21:347 | Construct CI on theta1 - R*theta2 not on R directly; allows asymmetric CIs |

### Odds ratio

| Method | Citation | Verdict |
|--------|----------|---------|
| Wald on log(OR) | textbook | Fast; suffers Hauck-Donner effect when cell counts small |
| Cornfield exact | Cornfield 1956 | Exact reference standard; often over-conservative |
| Mid-p Cornfield | Berry-Armitage 1995 *Statistician* 44:417 | Less conservative; slight under-coverage |
| Profile likelihood | Venzon-Moolgavkar 1988 *Appl Stat* 37:87 | **Transformation-invariant; no Hauck-Donner pathology**; R `MASS::confint.glm` default for `glm` |

**Hauck-Donner effect (1977 *JASA* 72:851; revived by Yee 2022 *JASA* 117:1763):** the Wald test statistic is *non-monotonic* in the parameter estimate near the boundary -- a large OR can produce a tiny Wald chi-square so the test fails to reject when it should. Use LR or profile-likelihood inference. `VGAM::hdeff()` detects this in fitted models.

## Number Needed to Treat (NNT) -- the Bender 2001 *CCT* 22:102 Way

**Goal:** Convert RD to clinically intuitive NNT with a CI that handles the singularity at RD = 0.

**Approach:** Compute RD CI first, then transform; when RD CI crosses zero, report NNTB(lower) -> infinity -> NNTH(upper).

```python
import numpy as np

def nnt_with_ci(treated_events, treated_n, control_events, control_n, alpha=0.05):
    """Bender 2001 *CCT* 22:102 method; transforms RD CI to NNT CI handling singularity."""
    p_t = treated_events / treated_n
    p_c = control_events / control_n
    rd = p_c - p_t   # positive = treatment helps
    se_rd = np.sqrt(p_t * (1 - p_t) / treated_n + p_c * (1 - p_c) / control_n)
    z = 1.96 if alpha == 0.05 else None
    rd_ci = (rd - z * se_rd, rd + z * se_rd)

    nnt = 1 / rd if rd != 0 else float('inf')
    if rd_ci[0] > 0 and rd_ci[1] > 0:
        return f'NNTB {1/rd_ci[1]:.0f} to {1/rd_ci[0]:.0f}'
    elif rd_ci[0] < 0 and rd_ci[1] < 0:
        return f'NNTH {1/-rd_ci[1]:.0f} to {1/-rd_ci[0]:.0f}'
    else:
        # CI crosses zero -- Bender convention
        if rd > 0:
            return f'NNTB {nnt:.0f} (NNTB {1/rd_ci[1]:.0f} to inf to NNTH {1/-rd_ci[0]:.0f})'
        return f'NNTH {-nnt:.0f} (similar split)'
```

**Bender 2001 *Controlled Clinical Trials* 22:102-110** resolves the discontinuity at RD = 0: NNT has a singularity there, so a CI that crosses zero produces a disjoint NNT CI ("NNTB(some) -> infinity -> NNTH(some)"). The Cochrane/BMJ convention is to report exactly this -- the infinity in the middle signals non-significance and is more honest than truncating to one side.

**NNT from OR + baseline risk** (when only the OR is published):

```python
def nnt_from_or(odds_ratio, baseline_risk):
    baseline_odds = baseline_risk / (1 - baseline_risk)
    treatment_odds = baseline_odds * odds_ratio
    treatment_risk = treatment_odds / (1 + treatment_odds)
    arr = abs(baseline_risk - treatment_risk)
    return 1 / arr if arr > 0 else float('inf')
```

| OR | Baseline 5% | Baseline 20% | Baseline 50% |
|----|-------------|--------------|--------------|
| 0.5 | NNT=42 | NNT=12 | NNT=6 |
| 0.7 | NNT=70 | NNT=20 | NNT=12 |

**Always report baseline risk alongside NNT** -- the same OR produces dramatically different NNTs.

## Marginal RD via G-Computation -- The FDA 2023 Recipe

**Goal:** Compute the marginal RD (FDA 2023 primary estimand for binary endpoints) from a fitted logistic regression -- without requiring the analyst to switch model class.

**Approach:** Standardise over the observed covariate distribution. Predict per-subject probability under treatment=1 AND under treatment=0; take the mean of each; the difference is the marginal RD. SE via influence function or bootstrap.

```python
import statsmodels.formula.api as smf
import numpy as np

# Step 1: Fit logistic with covariates
fit = smf.logit('y ~ z + x1 + x2 + x3', data=df).fit()

# Step 2: Predict under each treatment regime (counterfactual prediction)
df_z1 = df.assign(z=1)
df_z0 = df.assign(z=0)
p_z1 = fit.predict(df_z1)
p_z0 = fit.predict(df_z0)

# Step 3: Marginal estimates
marg_p1 = p_z1.mean()
marg_p0 = p_z0.mean()
marg_rd = marg_p1 - marg_p0          # FDA 2023 primary
marg_rr = marg_p1 / marg_p0
marg_or = (marg_p1 / (1 - marg_p1)) / (marg_p0 / (1 - marg_p0))

# Step 4: SE via bootstrap (or analytical influence function in R `marginaleffects`)
n_boot = 1000
boot_rds = np.zeros(n_boot)
rng = np.random.default_rng(42)
for b in range(n_boot):
    idx = rng.choice(len(df), size=len(df), replace=True)
    fit_b = smf.logit('y ~ z + x1 + x2 + x3', data=df.iloc[idx]).fit(disp=0)
    p1_b = fit_b.predict(df_z1.iloc[idx]).mean()
    p0_b = fit_b.predict(df_z0.iloc[idx]).mean()
    boot_rds[b] = p1_b - p0_b

se_marg_rd = boot_rds.std(ddof=1)
ci_marg_rd = (marg_rd - 1.96*se_marg_rd, marg_rd + 1.96*se_marg_rd)
```

**R equivalent (preferred for confirmatory):**

```r
library(marginaleffects)
fit <- glm(y ~ z + x1 + x2 + x3, family = binomial, data = df)
avg_comparisons(fit, variables = 'z', vcov = 'HC3')
# Returns marginal RD with HC3 sandwich SE; CIs; deltamethod_inference
```

**Tsiatis et al 2008 robustness guarantee:** under randomisation Z ⊥ X, the g-computation marginal estimator is consistent for marginal ATE EVEN IF the outcome model is misspecified. Efficiency depends on model quality; consistency does not. This is why FDA 2023 accepts marginal RD via g-computation without requiring proof of correct logistic mean structure.

## Marginal vs Conditional Effects -- The Core ICH E9(R1) Question

For logistic regression, the maximum-likelihood coefficient on Z in `glm(Y ~ Z + X, family=binomial)` is the **conditional log odds ratio** -- the OR comparing Z=1 to Z=0 *holding X fixed*. Because the OR is non-collapsible, this is *different* from the marginal log OR even under perfect randomisation. Lin & Wei 1989 first proved the ANCOVA-style estimator from a nonlinear model is generally inconsistent for the marginal ATE; this is the result FDA 2023 implicitly relies on.

```python
# Conditional OR (the default; biased toward null vs marginal OR)
import statsmodels.formula.api as smf
import numpy as np

model = smf.logit('y ~ z + x1 + x2', data=df).fit()
cond_or = np.exp(model.params['z'])

# Marginal RD via g-computation / standardisation:
df_z1, df_z0 = df.assign(z=1), df.assign(z=0)
p_z1 = model.predict(df_z1).mean()
p_z0 = model.predict(df_z0).mean()
marg_rd = p_z1 - p_z0
marg_rr = p_z1 / p_z0
marg_or = (p_z1 / (1-p_z1)) / (p_z0 / (1-p_z0))

# Variance via influence function or bootstrap (RobinCar / marginaleffects)
# R: marginaleffects::avg_comparisons(model, variables='z', vcov='HC3')
```

**The reporting standard post-FDA-2023:** marginal estimand as primary (RD with HC3 sandwich SE), conditional as supportive. The two are different parameters, not different estimates of the same parameter. See clinical-biostatistics/logistic-regression for full g-computation machinery and Tsiatis et al 2008 *Stat Med* 27:4658.

### Non-collapsibility -- the standard explanation

A conditional (adjusted) OR typically differs from the marginal (unadjusted) OR even WITHOUT confounding. Conditioning on a prognostic covariate (one that predicts outcome but is unrelated to treatment) changes the OR by a mathematical property -- not bias. In most settings, the conditional OR is *further from the null* than the marginal OR. **This is what postdocs argue about:** Frank Harrell argues conditional effects are more clinically meaningful (they apply to individuals); FDA via Permutt argues marginal effects are what regulatory submissions need (the population-level decision). FDA 2023 settles it for regulatory reporting: marginal primary, conditional supportive.

**Risk ratios and risk differences are collapsible.** Marginal and conditional RR/RD are equal absent confounding. This is why FDA 2023 favours RD as the primary estimand for binary outcomes -- the conditional vs marginal distinction disappears.

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| Conditional OR (from adjusted logistic) > marginal OR (from g-computation) | Non-collapsibility (Permutt 2020); adjusted ORs further from null than marginal under most covariate structures | Report marginal RD as primary per FDA 2023; conditional OR as supportive with explicit parameter label |
| Wald CI excludes null but Wilson/Newcombe-Wilson CI overlaps null | Wald has poor coverage near 0 and 1 (Brown-Cai-DasGupta 2001) | Wilson for single proportion; Newcombe-Wilson hybrid or Miettinen-Nurminen for differences; cite as regulatory standard |
| Crude OR vs stratified MH-OR differ substantially | Confounding by stratification factor OR non-collapsibility OR effect modification | If RCT with stratified randomisation, use stratified analysis (Kahan-Morris 2012); if observational, investigate confounding via stratum-specific OR forest plot |
| OR vs RR direction agree but magnitudes differ greatly | Outcome prevalence > 10% (OR overstates RR at high baseline risk) | Use modified Poisson (Zou 2004) for direct RR estimation; report both for transparency |
| NNT from one trial vs another differs despite similar OR | Baseline risk differs across trials (NNT = 1/ARR depends on baseline) | Always report NNT alongside baseline risk; cite Bender 2001 *CCT* 22:102 for NNTB-inf-NNTH convention |
| Profile-likelihood CI differs from Wald CI in small-cell scenario | Hauck-Donner effect (Wald non-monotone near boundary; Yee 2022) | Use profile likelihood (R `MASS::confint.glm`); detect via `VGAM::hdeff()` |
| RMST difference and HR give different conclusions on treatment benefit | HR is time-averaged log-HR under PH violation; RMST captures cumulative benefit | Under PH violation, RMST is the more interpretable summary; see clinical-biostatistics/survival-analysis |

## Modified Poisson Regression for Common Outcomes

When prevalence > 10% and the policy quantity is RR (not OR), modified Poisson with sandwich SE directly estimates the RR (Zou 2004 *AJE* 159:702):

```python
import statsmodels.api as sm
import numpy as np

poisson_model = sm.GLM(df['outcome'],
                       sm.add_constant(df[['treatment', 'age']]),
                       family=sm.families.Poisson()).fit(cov_type='HC1')
rr_estimates = np.exp(poisson_model.params)
rr_ci = np.exp(poisson_model.conf_int())
```

`cov_type='HC1'` (Huber-White sandwich) corrects the over-dispersion that Poisson assumes; without it, SEs are wrong because binary data are NOT Poisson. Use HC3 (jackknife approximation) for n < 250 (Long-Ervin 2000 *Am Stat* 54:217). HC1 is the Stata default; HC3 is the R `sandwich` package default -- the difference can flip noninferiority p-values for n < 200.

## Forest Plots

```python
import matplotlib.pyplot as plt
import numpy as np

def forest_plot(labels, effects, lower_cis, upper_cis, ref_line=1.0, figsize=(8, 6), log_scale=True):
    fig, ax = plt.subplots(figsize=figsize)
    y_pos = range(len(labels))
    ax.errorbar(effects, y_pos,
                xerr=[np.array(effects) - np.array(lower_cis),
                      np.array(upper_cis) - np.array(effects)],
                fmt='D', color='black', capsize=3, markersize=5)
    ax.axvline(x=ref_line, color='gray', linestyle='--', linewidth=0.8)
    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(labels)
    ax.set_xlabel('Effect estimate (95% CI)')
    if log_scale:
        ax.set_xscale('log')
    plt.tight_layout()
    return fig
```

Log scale ensures reciprocal effects (OR 0.5 and OR 2.0) appear equidistant from the null. For RD or RMST difference, use a linear scale.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Prevalence > 10% -> prefer RR over OR | Zhang-Yu 1998 *JAMA* 280:1690 | OR overstates RR materially; modified Poisson or log-binomial recovers RR directly |
| n > 100 per arm for Wald CI on log-OR | Brown-Cai-DasGupta 2001; Newcombe 1998a | Below this, Wald coverage is poor; profile likelihood or score-based preferred |
| HC1 (Stata default) vs HC3 (R default) for small n | Long-Ervin 2000 *Am Stat* 54:217 | HC3 (jackknife approximation) recommended for n <=250; difference can flip NI p-values |
| Marginal estimand for primary regulatory report | FDA 2023 Final Guidance | Marginal RD/RR/OR via g-computation; conditional OR is a different parameter (Permutt 2020) |
| Miettinen-Nurminen CI for RR/RD | Miettinen-Nurminen 1985; FDA/EMA NI margin practice | Regulatory standard; consistent with Pearson chi-square |
| Bender NNT convention | Bender 2001 *Controlled Clinical Trials* 22:102-110; Cochrane/BMJ style | When RD CI crosses zero, NNT CI is disjoint; report NNTB(lower) -> inf -> NNTH(upper) |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Adjusted OR much larger than unadjusted | Non-collapsibility, not confounding | Cite Permutt 2020; report both with explicit estimand labels (conditional vs marginal) |
| NNT reported as "NNT=25, CI 12 to -200" | Sign-confused output when RD CI crosses zero | Use Bender 2001 *CCT* 22:102 convention: "NNTB 25 (NNTB 12 to inf to NNTH 200)" |
| OR reported without baseline risk for clinical translation | Common in published papers | Always report event rates per arm alongside OR; provide NNT at observed baseline |
| `Table2x2` returns reciprocal OR | Column ordering puts outcome=0 first | Reorder: `crosstab[[1, 0]]` |
| Poisson SE much smaller than expected | Forgot `cov_type='HC1'` or 'HC3' | Always specify sandwich SE for modified Poisson; without it, SEs are wrong |
| Wald log-RR CI extends below 0 or implausibly | Katz log RR fails for small p | Switch to Miettinen-Nurminen or Koopman score; R `ratesci::scoreci(contrast='RR')` |
| Hauck-Donner: huge OR with non-significant Wald p | Wald non-monotonic near boundary | Use profile likelihood: R `MASS::confint.glm`; detect with `VGAM::hdeff()` |
| "Adjusted" model gives smaller effect than expected | Mediator adjustment | Check causal DAG; adjusting for mediators attenuates effect toward null |

## Anticipated Reviewer Pushback

| Pushback | Response |
|----------|----------|
| "Is the OR collapsible across the included covariates?" | No -- OR is non-collapsible. Cite Permutt 2020. Report marginal RD via g-computation as primary per FDA 2023. |
| "Why MN over Wald CI?" | Wald has poor coverage for sparse 2x2; MN matches Pearson chi-square test and is the FDA/EMA standard for NI margins. |
| "Hauck-Donner pathology check?" | For small cell counts: switch to profile likelihood CI (R `MASS::confint.glm`); cite Yee 2022 *JASA* 117:1763. |
| "Why HC3 over HC1?" | n <=250 favours HC3 per Long-Ervin 2000; HC1 (Stata) and HC3 (R sandwich) can disagree at small n. |
| "Where is the marginal effect?" | Per FDA 2023, marginal RD/RR via g-computation is the primary estimand for binary outcomes; conditional is supportive. |
| "Why NNT in 'NNTB-infinity-NNTH' notation?" | Bender 2001 *CCT* 22:102 convention; standard in BMJ/Lancet/Cochrane. Disjoint CI honestly conveys non-significance. |
| "Adjustment for stratification factors?" | Strata included in modified Poisson or in MN stratified CI via `ratesci::scoreci(stratified=TRUE)`; ignoring inflates Type-I (Kahan-Morris 2012). |
| "What about effect modification across subgroups?" | Pooled estimate reported as primary; stratum-specific ORs in forest plot; Breslow-Day test for homogeneity; if heterogeneous, do not pool -- see clinical-biostatistics/subgroup-analysis. |
| "Post-hoc subgroup OR was significant -- can a claim be made?" | Not credible without pre-specification per EMA 2019 / CONSORT 2025; frame as hypothesis-generating; replicate in independent cohort before claiming. |

## References

- Agresti A, Caffo B. 2000. Simple and effective confidence intervals for proportions and differences of proportions. *Am Stat* 54:280-288.
- Altman DG. 1998. Confidence intervals for the number needed to treat. *BMJ* 317:1309.
- Bender R. 2001. Calculating confidence intervals for the number needed to treat. *Controlled Clinical Trials* 22:102-110.
- Brown LD, Cai TT, DasGupta A. 2001. Interval estimation for a binomial proportion. *Stat Sci* 16:101-117.
- Chan ISF, Zhang Z. 1999. Test-based exact confidence intervals for the difference of two binomial proportions. *Biometrics* 55:1202-1209.
- Cornfield J. 1956. A statistical problem arising from retrospective studies. *3rd Berkeley Symp* 4:135-148.
- Donner A, Zou GY. 2012. Closed-form confidence intervals for ratios involving population means via MOVER-R. *Stat Methods Med Res* 21:347-359.
- FDA. 2023. Adjusting for Covariates in Randomized Clinical Trials. Final Guidance, May 2023.
- Hauck WW, Donner A. 1977. Wald's test as applied to hypotheses in logit analysis. *JASA* 72:851-853.
- Katz D, Baptista J, Azen SP, Pike MC. 1978. Obtaining confidence intervals for the risk ratio in cohort studies. *Biometrics* 34:469-474.
- Koopman PAR. 1984. Confidence intervals for the ratio of two binomial proportions. *Biometrics* 40:513-517.
- Long JS, Ervin LH. 2000. Using heteroscedasticity consistent standard errors in the linear regression model. *Am Stat* 54:217-224.
- Miettinen O, Nurminen M. 1985. Comparative analysis of two rates. *Stat Med* 4:213-226.
- Newcombe RG. 1998a. Interval estimation for the difference between independent proportions: comparison of eleven methods. *Stat Med* 17:873-890.
- Permutt T. 2020. Do covariates change the estimand? *Stat Biopharm Res* 12:45-53.
- Royston P, Parmar MKB. 2013. Restricted mean survival time. *BMC Med Res Methodol* 13:152.
- Tsiatis AA, Davidian M, Zhang M, Lu X. 2008. Covariate adjustment for two-sample treatment comparisons in randomized clinical trials. *Stat Med* 27:4658-4677.
- Venzon DJ, Moolgavkar SH. 1988. A method for computing profile-likelihood-based confidence intervals. *Appl Stat* 37:87-94.
- Yee TW. 2022. On the Hauck-Donner effect in Wald tests: detection, tipping points, and parameter space characterization. *JASA* 117:1763-1774.
- Zou G. 2004. A modified Poisson regression approach to prospective studies with binary data. *AJE* 159:702-706.

## Related Skills

- clinical-biostatistics/categorical-tests - 2x2 testing infrastructure that produces ORs/RRs/RDs
- clinical-biostatistics/logistic-regression - Adjusted ORs, g-computation, modified Poisson for marginal RR
- clinical-biostatistics/subgroup-analysis - Forest plots and stratified effect estimates
- clinical-biostatistics/survival-analysis - HR, RMST, hazard-free alternatives for time-to-event
- clinical-biostatistics/trial-reporting - CONSORT 2025 and ICH E9(R1) effect reporting
- clinical-biostatistics/missing-data-sensitivity - Effect measure CIs under MI pooling (Rubin's rules)

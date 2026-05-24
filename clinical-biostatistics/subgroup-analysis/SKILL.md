---
name: bio-clinical-biostatistics-subgroup-analysis
description: Performs subgroup and heterogeneous treatment effect (HTE) analyses for clinical trials. Covers Mantel-Haenszel pooling, Breslow-Day, interaction tests in regression, RERI for additive interaction, modern data-adaptive HTE methods (STEPP, SIDES, causal forests, X/R-learners), Bayesian shrinkage (Dixon-Simon, MAP, EXNEX), graphical multiplicity (Bretz-Maurer), and credibility frameworks (Sun BMJ, EMA 2019). Use when analyzing treatment effects across patient subgroups for regulatory submissions or precision-medicine claims.
tool_type: python
primary_tool: statsmodels
---

## Version Compatibility

Reference examples tested with: statsmodels 0.14+, scipy 1.12+, numpy 1.26+, pandas 2.1+, matplotlib 3.8+, scikit-learn 1.4+. R packages cited: grf, policytree, causalToolbox, personalized, SIDES, stepp, gMCP, partykit, RBesT, brms.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name`

If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# Subgroup Analysis and Heterogeneous Treatment Effects

**"Analyze treatment effects across subgroups"** -> Test whether treatment effects differ across pre-specified or data-discovered subgroups using interaction tests, stratified estimators, modern data-adaptive HTE methods, or Bayesian shrinkage -- with explicit declaration of confirmatory vs exploratory intent and credibility assessment.

## The Senn Foundation -- Why Most Subgroup Claims Are Wrong

**Senn 2018 *Nature* 563:619-621 (and *Statistical Issues in Drug Development* Ch. 9, 14):** observed between-patient response variation is NOT evidence of patient-level HTE. It conflates within-patient noise, period effects, regression-to-the-mean, and measurement error with true individual heterogeneity. Senn-Rolfe-Julious 2011 *SMMR* 20:657 documents that variance-component decomposition of replicate-crossover trials repeatedly fails to find subject-by-treatment interaction even where reviewers were certain one must exist.

**Brookes et al 2004 *J Clin Epidemiol* 57:229 — the 4x penalty:** detecting a treatment-by-subgroup interaction requires approximately 4x the sample size needed to detect the main treatment effect of similar magnitude. A trial powered to detect OR=0.6 overall cannot reliably detect subgroup differences of similar magnitude. Non-significant interaction tests are usually underpowered, not null.

**Senn's aphorism (paraphrased from SIDD):** "a trial can have subgroup analyses or proper power, not both." A single trial cannot simultaneously be powered for a primary effect AND for credible subgroup discovery; pretending otherwise misrepresents posterior uncertainty.

## Algorithmic Taxonomy

| Method | What it answers | Inference | Strength | Fails when |
|--------|-----------------|-----------|----------|------------|
| Mantel-Haenszel / CMH | Common OR across pre-defined strata | Asymptotic | Preserves stratification factor from randomisation | Stratum ORs reverse direction (Simpson) |
| Breslow-Day | Homogeneity of stratum ORs | Asymptotic chi-square | Test of effect modification | Underpowered with few/sparse strata; non-significance NOT proof of homogeneity |
| Gail-Simon 1985 | **Qualitative** interaction (sign reversal) | Likelihood ratio | Distinguishes quantitative from qualitative | Original LR is liberal at small n; use exact critical values (Pan-Wolfe 1997) |
| Logistic regression interaction | Per-subgroup conditional OR | Wald, LR | Most efficient single-model approach | Conditional ORs subject to non-collapsibility; over-fits in multi-way |
| RERI | Additive interaction on multiplicative scale | Delta-method or bootstrap CI | Captures public-health-relevant scale | Delta-method poor near boundary; nonlinear function of three ORs |
| STEPP (Bonetti-Gelber 2000) | Continuous covariate subgroups via overlapping windows | Permutation supremum test | Avoids dichotomisation of biomarkers | Window-size choice affects results; correlated estimates require permutation inference |
| SIDES / SIDEScreen (Lipkovich 2011) | Data-discovered subgroups via recursive partitioning | Resampling-adjusted base-vs-complement p | Multiplicity correctly absorbed | Tuning skeleton parameters affects FWER calibration |
| QUINT (Dusseldorp 2014) | Qualitative-interaction trees | Bootstrap stability | Directly tests crossover (A-better, B-better, equal) | Only two-arm continuous/binary; survival extensions ad hoc |
| Virtual Twins (Foster 2011) | Per-subject CATE via twin RF predictions | Bootstrap | Decouples nuisance from interpretable subgroup | Biased when RF underestimates effect heterogeneity in either arm |
| Causal forests (Athey-Wager 2019) | Pointwise CATE with honest splits | Influence-function CI | Asymptotic Gaussianity; doubly robust via AIPW | Finite-sample CI validity at trial-scale n debated (Rehill 2025) |
| Meta-learners X/R-learner (Künzel 2019) | Marginal/conditional CATE | Cross-fit influence function | X-learner dominates T-learner when arms unbalanced | Needs propensity for X-learner; R-learner requires nuisance n^(1/4) rate |
| MOB (Zeileis 2008) | Parameter-instability trees | M-fluctuation test | Theoretically clean; invariant to monotone transforms | Worse out-of-sample CATE than causal forest at large p |
| Bayesian shrinkage (Dixon-Simon 1991) | Posterior subgroup effects shrunken to overall | Posterior intervals from MCMC | Honest about prior expectation of no qualitative interaction | Prior choice on tau drives results; Dane et al 2019 white paper warns against for signal generation |
| EXNEX (Neuenschwander 2016) | Mixture of exchangeable + per-basket non-exchangeable | Posterior | Avoids HM "catastrophic borrowing" when one basket truly different | Weight choice (often 0.5/0.5 default) affects borrowing strength |

**Postdoc reading list:**

- Wang et al 2007 *NEJM* 357:2189 ("Reporting of subgroup analyses in clinical trials") — the canonical NEJM-mandated practice
- Sun et al 2010/2012 *BMJ* 340:c117 and 344:e1553 — 11 credibility criteria
- Dane, Spencer, Rosenkranz, Lipkovich, Parke 2019 *Pharm Stat* 18:126 with Hemmings-Koch commentary at 18:140 — EFSPI white paper + critique
- EMA 2019 Guideline on subgroups (CHMP/EWP/117211/2010 effective Aug 2019) — distinguishes "assessment subgroups" (regulatory, pre-specified) from "discovery subgroups" (exploratory)
- Athey & Wager 2019 *Observational Studies* 5:37; Athey, Tibshirani, Wager 2019 *Ann Stat* 47:1148 — causal forests
- Rehill 2025 *Int Stat Rev* — applied causal forest audit; ~70% of papers misreport tuning and skip honest-splitting validation

## Decision Tree by Scenario

| Scenario | Recommended approach | Why |
|----------|----------------------|-----|
| Pre-specified subgroup in confirmatory trial | Interaction term in single model + graphical-procedure multiplicity (gMCP) | EMA 2019 "assessment subgroup"; interaction test + Bonferroni-graph alpha control |
| Stratified randomisation by site/region | CMH or logistic with strata as covariates; stratum-specific OR reported in forest plot | Kahan-Morris 2012; ignoring strata inflates Type-I; strata-specific ORs for transparency |
| 5-15 pre-specified subgroups, regulatory submission | Forest plot + interaction tests + Holm/graphical FWER correction | Default regulatory presentation; multiplicity adjustment expected |
| Continuous biomarker subgroup (e.g., HbA1c, biomarker score) | STEPP (sliding-window plot + permutation supremum test) | Avoids arbitrary dichotomisation; cite Bonetti-Gelber 2004 |
| Suspected qualitative interaction (treatment helps some, harms others) | Gail-Simon 1985 LR test with Pan-Wolfe 1997 exact critical values | Distinguishes quantitative from qualitative; critical for label restriction |
| Data-discovery of HTE subgroup | SIDES/SIDEScreen with permutation FWER + replication mandate | Lipkovich 2011; signal-discovery not signal-confirmation |
| Continuous CATE estimation with many covariates | Causal forest (R `grf::causal_forest`) + RATE test (Yadlowsky 2025) | Modern HTE; honest splitting; influence-function CIs; RATE tests whether ranking is predictive vs prognostic |
| Basket trial across rare-disease strata | EXNEX or robust MAP with `RBesT` | Borrows across strata while permitting one to detach if truly different |
| Subgroup signal needing replication planning | Bayesian shrinkage for adjusted estimates; Sun 2010 winner's curse | Selected subgroups have inflated effect by selection bias |
| Pediatric extrapolation borrowing from adult data | Power prior (gamma in 0.3-0.6) per FDA Bayesian draft Jan 2026 | Partial borrowing with discount; standard regulatory approach |

## Mantel-Haenszel and Stratified Analysis

**Goal:** Estimate a pooled treatment effect across strata while testing homogeneity.

**Approach:** Construct per-stratum 2x2 tables; CMH for pooled OR and null test; Breslow-Day for homogeneity (with caveats).

```python
from statsmodels.stats.contingency_tables import StratifiedTable
import pandas as pd
import numpy as np

tables = []
for stratum in df['site'].unique():
    sub = df[df['site'] == stratum]
    t = pd.crosstab(sub['treatment'], sub['outcome']).values
    if t.shape == (2, 2) and t.min() > 0:
        tables.append(t)

st = StratifiedTable(tables)
print('MH pooled OR:', st.oddsratio_pooled)
print('95% CI:', st.oddsratio_pooled_confint())
print('CMH H0: common OR=1, p =', st.test_null_odds().pvalue)
print('Breslow-Day H0: equal ORs, p =', st.test_equal_odds().pvalue)
```

**Breslow-Day power trap:** with k=3 strata and modest heterogeneity, Breslow-Day power can be <40%. Non-significance does NOT prove homogeneity — it just means the null cannot be rejected. Always supplement with a forest plot AND an LR interaction test from logistic regression.

## Interaction Terms in Regression -- The Correct Way

**Single model with interaction is the regulatory standard** — comparing p-values from separate per-subgroup models is statistically invalid (separate models have different power, and the p-value differences confound effect size with sample size).

```python
import statsmodels.formula.api as smf
import numpy as np

# Single model with interaction (CORRECT)
model = smf.logit(
    'outcome ~ C(treatment, Treatment(reference="Placebo")) * C(age_group)',
    data=df
).fit()
# Interaction coefficient tests effect modification
# LR test: compare to additive model:
additive = smf.logit('outcome ~ C(treatment, Treatment(reference="Placebo")) + C(age_group)', data=df).fit()
lr_p = 1 - chi2.cdf(2 * (model.llf - additive.llf), model.df_model - additive.df_model)

# Extract subgroup-specific ORs for reporting:
for group in df['age_group'].unique():
    sub_model = smf.logit('outcome ~ C(treatment, Treatment(reference="Placebo"))',
                          data=df[df['age_group'] == group]).fit()
    or_val = np.exp(sub_model.params.iloc[1])
    ci = np.exp(sub_model.conf_int().iloc[1])
    print(f'{group}: OR={or_val:.3f} ({ci[0]:.3f}-{ci[1]:.3f})')
```

## RERI for Additive Interaction

```python
# Fit interaction model: outcome ~ treatment + subgroup_indicator + treatment:subgroup_indicator
# OR_11 = OR for treated in subgroup (vs untreated not in subgroup)
# OR_10 = OR for treated not in subgroup
# OR_01 = OR for untreated in subgroup
reri = or_11 - or_10 - or_01 + 1
# RERI > 0 = synergism (combined effect > sum of individual)
# RERI = 0 = no additive interaction
# RERI < 0 = antagonism
# CI requires delta method or bootstrap (nonlinear function of three ORs)
```

**Multiplicative vs additive interaction:** logistic regression tests multiplicative interaction (ratio of ORs). Null multiplicative interaction does NOT imply null additive interaction. For public health decisions, additive interaction is often more relevant.

## Modern Data-Adaptive HTE Methods

### STEPP (Subpopulation Treatment Effect Pattern Plot)

```python
# R recommended; Python equivalents are emerging
# library(stepp)
# subset_obj <- new('stwin', type='sliding', r1=100, r2=200)  # window sizes
# step_obj <- stepp(eff='binary', cov='biomarker', trt='treatment',
#                   resp='outcome', subset=subset_obj, ...)
# plot(step_obj); test_pattern(step_obj, nperm=1000)
```

Bonetti-Gelber 2000/2004; window-size choice (sliding vs tail-oriented) affects results. **Naive simultaneous CIs are wrong** — estimates are correlated across windows; use permutation-based supremum tests of pattern flatness (Yip et al 2016 *Stat Med* 35:5437 documents 30% disagreement between sliding and tail-oriented variants).

### SIDES / SIDEScreen

```python
# R: library(SIDES) or library(rsides)
# Recursive partitioning with differential-effect splitting + permutation-adjusted subgroup p
```

Lipkovich-Dmitrienko-Denne-Enas 2011 *Stat Med* 30:2601; SIDEScreen (Lipkovich-Dmitrienko 2014 *J Biopharm Stat* 24:130) adds variable-importance prefilter + base-vs-complement test. **The inferential complement test correctly absorbs multiplicity of considered splits** — Bonferroni alternatives over-correct.

### Causal forests (Athey-Wager)

```python
# R recommended; econml is the Python equivalent
# library(grf)
# cf <- causal_forest(X, Y, W, num.trees=2000)
# tau.hat <- predict(cf, X)$predictions
# test_calibration(cf)  # omnibus calibration test
```

**Honest splitting** (one sub-sample for splits, another for leaf estimates) yields asymptotic Gaussianity and pointwise CIs (Athey-Tibshirani-Wager 2019 *Ann Stat* 47:1148). Doubly-robust variants use AIPW pseudo-outcomes.

**Diagnostic discipline (Rehill 2025 audit):** ~70% of applied causal-forest papers misreport tuning, skip honest-splitting validation, or omit the `test_calibration` check. Required diagnostic steps:

1. Honest splitting enabled (`honesty=TRUE`)
2. `test_calibration(cf)` — regresses actual treatment effects on out-of-bag CATE predictions; significant positive slope = CATE has signal
3. Variable importance via permutation
4. RATE/AUTOC test (Yadlowsky 2025 *JASA*) — single p-value omnibus test for whether CATE ranking has predictive (not just prognostic) value

### Meta-learners (Künzel-Sekhon-Bickel-Yu 2019 *PNAS* 116:4156)

- **S-learner:** single outcome model on (Z, X); take difference. Biased when treatment effect smaller than baseline effect.
- **T-learner:** separate models per arm; take difference. Standard baseline.
- **X-learner:** T-learner + cross-imputation + propensity weighting. Dominates T-learner when arms unbalanced (1:k, k>=3) or true CATE is smoother than response.
- **R-learner (Nie-Wager 2021 *Biometrika* 108:299):** orthogonal Robinson-style residualisation, debiased ML. Asymptotically dominates X-learner when nuisance models converge at n^(1/4) rate.

```python
# Python: econml
from econml.dml import CausalForestDML, LinearDML
from econml.metalearners import XLearner, TLearner
# Or: from causalml.inference.meta import LRSRegressor, XGBTRegressor

xl = XLearner(models=GradientBoostingRegressor(), propensity_model=LogisticRegression())
xl.fit(Y, T=W, X=X)
cate = xl.effect(X)
```

## Bayesian Shrinkage -- The Postdoc Argument

**Dixon-Simon 1991 *Biometrics* 47:871:** exchangeable mean-zero prior on treatment-by-subgroup interaction coefficients; subgroup posteriors shrink to overall trial effect, proportional to evidence of heterogeneity. Honest about prior expectation that no qualitative interaction will hold.

**Berry, Broglio, Groshen, Berry 2013 *Clin Trials* 10:720 (NOT JCO — common citation error):** hierarchical pooling across baskets calibrated by between-basket variance tau-squared. Inflated false-positives when tau-squared mis-specified (Freidlin-Korn 2013 *JNCI* 105:1532 critique).

**EXNEX (Neuenschwander, Wandel, Roychoudhury, Bailey 2016 *Pharm Stat* 15:123):** mixture of exchangeable (shared mean+variance) + non-exchangeable (per-basket independent) components, typically weighted 0.5/0.5. Avoids HM catastrophic borrowing when one basket truly different. Now near-standard in oncology basket trials.

**MAP priors (Schmidli et al 2014 *Biometrics* 70:1023):** meta-analytic-predictive prior. Fit random-effects meta-analysis of historical control arms, derive predictive distribution for new control arm, use as informative prior. Effective sample size from history typically 20-80% of new control arm. **Robust MAP** adds vague mixture component (weight 0.1-0.3) to guard against prior-data conflict. R `RBesT` is the canonical CRAN tool (Weber et al 2021 *J Stat Softw* 100:19).

**The Dane vs Hemmings argument:**

- Dane et al 2019 *Pharm Stat* 18:126 (EFSPI white paper): propose a standardised-effect benchmark in the neighbourhood of ~0.5σ for "noteworthy heterogeneity" (verify the exact figure against the white paper before quoting); tiered structure (key/important/exploratory); signal-vs-noise diagnostics.
- Hemmings & Koch 2019 *Pharm Stat* 18:140 commentary: **explicitly REJECT Bayesian shrinkage for signal generation** because "Bayesian shrinkage assumes treatment effects are consistent" and pre-emptively damps the very heterogeneity one is searching for. Shrinkage is endorsed only for *post-signal replication planning*.

**What postdocs argue about:** whether shrinkage is appropriate at the signal-generation stage; the prior on tau drives everything (Senn-style: tau ~ HalfNormal(0, 0.1); regulatory-tolerant: tau ~ HalfNormal(0, 0.5)).

## Multiplicity for Subgroup Analyses

**Bonferroni is rarely right** for subgroup tests because they are heavily positively correlated (same outcome, overlapping samples). Bonferroni assumes worst-case dependence; loses 30-50% power vs Hommel/resampling for ~10 demographic subgroups.

**Graphical procedures (Bretz-Maurer-Hommel)** — see clinical-biostatistics/multiplicity-graphical for full treatment. The graph for a typical subgroup analysis SAP:

- Primary endpoint at full alpha
- Alpha propagates to "key secondary" endpoints if primary rejects
- Alpha propagates to pre-specified subgroup interaction tests if primary rejects
- Discovery subgroups at most receive a small alpha slice (Dane et al recommend <=20%)

```python
# R: library(gMCP); graphView(graphMCP(...))  # interactive or programmatic
```

**Goeman, Hemerik, Solari 2021 *Ann Stat* 49:1218:** "only closed testing procedures are admissible for controlling FDP/FWER/k-FWER" — graphical, gatekeeping, Hommel, fallback are all closed tests in disguise. The graph IS the procedure.

## Pre-Specified vs Post-Hoc and EMA 2019

| Aspect | Pre-specified ("assessment") | Post-hoc ("discovery") |
|--------|------------------------------|------------------------|
| Timing | Before unblinding, in SAP | After seeing data |
| Credibility | High if biologically justified | Low; hypothesis-generating only |
| Regulatory weight per EMA 2019 | Can support labeling claims | Cannot support claims alone |
| Multiplicity adjustment | Required per SAP (graphical or Holm) | Required + heavy skepticism |
| Number expected | Few (5-15); pre-justified | Unlimited but unblindable |

**EMA 2019 Guideline on Investigation of Subgroups in Confirmatory Trials (CHMP/EWP/117211/2010, effective Aug 2019)** position: interaction tests alone are "neither necessary nor sufficient"; consistency is a holistic judgment; subgroup-specific licensing requires pre-specified evidence of differential effect PLUS biological rationale.

**Sun BMJ 2012 11 credibility criteria** (the canonical academic framework):

1. Covariate measured pre-randomisation
2. A priori hypothesis
3. Direction pre-specified
4. One of few tests with multiplicity adjustment
5. Within-study comparison (not between)
6. Interaction test significant
7. Statistically independent subgroup
8. Large effect and all subgroups reported
9. Consistent across related outcomes
10. Consistent across studies
11. Biological plausibility / indirect evidence

## Quantitative vs Qualitative Interaction

**Gail-Simon 1985 *Biometrics* 41:361:** LR test against the "all-same-direction" null. Distinguishes quantitative (magnitude varies, direction preserved) from qualitative (sign reverses) — only the latter carries decision-changing weight clinically.

**Pan-Wolfe 1997 *Biometrics* 53:1564, Li-Chan 2006 *Stat Med* 25:2099:** exact critical values for Gail-Simon; original normal-approximation is liberal at small n.

**Why qualitative interaction matters for regulators:** may warrant restricting indication to the benefiting subgroup (label carve-out).

## The Winner's Curse and Sun et al 2012

Sun et al 2012 *BMJ* systematic review documented that the median observed effect size in *significant* subgroups across 64 trials was approximately 2.4x larger than the trial-overall effect — textbook winner's curse signature.

**Mechanism:** selecting a subgroup conditional on its interaction p being small selects for upward sampling fluctuations; the posterior MLE conditional on selection is biased upward by ~SE × inverse Mills ratio at the selection threshold.

**Mitigations:**

- Bayesian shrinkage (Dixon-Simon, RBesT) — posterior subgroup means shrink to overall
- Cross-validation (Athey-Wager honest splitting)
- Yadlowsky RATE 2025 *JASA* — single-p omnibus test for whether CATE ranking has predictive vs prognostic value

## Reconciliation: When Methods Disagree

| Pattern | Likely cause | Action |
|---------|--------------|--------|
| Causal forest CATE shows HTE; interaction test in logistic non-significant | Causal forest detects nonlinear/multivariate HTE; interaction test only linear bivariate | Causal forest with Yadlowsky RATE 2025 omnibus test is more sensitive; cite as discovery, not confirmatory; replicate in independent sample |
| STEPP pattern significant at one window choice; flat at another | Window-size choice (sliding vs tail-oriented) affects results (Yip 2016 documents 30% disagreement) | Pre-specify window choice in SAP; report sensitivity over choices; use permutation supremum test |
| Bayesian shrinkage (Dixon-Simon) yields null subgroup effect; frequentist sees signal | Shrinkage pulls subgroup posteriors toward overall trial effect | Bayesian shrinkage appropriate for replication PLANNING (Hemmings-Koch 2019), not signal generation; report frequentist as primary in discovery |
| MH stratified analysis gives different pooled OR than logistic with strata as covariates | MH conditions on stratum; logistic does not (assumes additivity of strata effects) | Both are valid under different assumptions; logistic preferred when stratum-treatment interaction tested formally |
| EXNEX basket trial gives different per-basket estimate under default 0.5/0.5 vs 0.3/0.7 mixture | Mixture weight controls borrowing strength (Neuenschwander 2016) | Sensitivity over weights 0.1-0.9; primary at 0.5/0.5; cite range |
| Gail-Simon qualitative interaction test rejects; LR interaction test does not | Qualitative test detects sign reversal; LR test detects magnitude | Both tests are answering different questions; qualitative interaction is regulatorily significant (label restriction) |
| Subgroup signal from data-adaptive HTE method (SIDES, causal forest) without replication | Winner's curse (Sun 2010); ~2.4x effect inflation in selected subgroups | Apply Bayesian shrinkage OR honest cross-validation; report as discovery only; require Phase 3 replication |
| Causal forest test_calibration passes; RATE/AUTOC fails | test_calibration measures any signal (prognostic OR predictive); RATE measures predictive value | RATE 2025 is the modern omnibus; should replace test_calibration as primary HTE check |

## Per-Method Failure Modes

### CMH masking Simpson's paradox

- **Trigger:** Stratum-specific ORs reverse direction; pooled MH OR appears null.
- **Mechanism:** CMH pools weighted log-ORs; opposite-sign equal-magnitude cancel.
- **Symptom:** Breslow-Day p < 0.05 with visible stratum sign reversal.
- **Fix:** Report stratum-specific as primary; switch to logistic with interaction term; pooled OR invalid.

### Breslow-Day false reassurance

- **Trigger:** Few strata (k<5) or sparse strata.
- **Mechanism:** Low power; chi-square with k-1 df.
- **Symptom:** Breslow-Day p>0.5 with forest plot showing visible heterogeneity.
- **Fix:** Forest plot + LR interaction test from logistic regression.

### STEPP correlated estimates with naive CI

- **Trigger:** STEPP plot with naive pointwise CIs at each window.
- **Mechanism:** Overlapping windows share patients; estimates are correlated.
- **Symptom:** Apparent significance at individual windows that disappears under permutation.
- **Fix:** Permutation-based supremum test of pattern flatness; report supremum p, not pointwise.

### Causal forest without honest splitting

- **Trigger:** Default settings without `honesty=TRUE`.
- **Mechanism:** In-sample split selection biases leaf estimates.
- **Symptom:** Over-fit CATE; failed `test_calibration`.
- **Fix:** Always `honesty=TRUE`; tune via `tune.parameters="all"`; cite Rehill 2025 audit.

### Bayesian shrinkage damping real heterogeneity

- **Trigger:** Hierarchical model fit during signal discovery (not replication planning).
- **Mechanism:** Prior on tau forces subgroup effects toward overall.
- **Symptom:** Subgroup detected by causal forest gets shrunken to null in shrinkage analysis.
- **Fix:** Hemmings-Koch position — shrinkage for replication planning, NOT signal generation; cite Dane et al 2019 + critique.

### EXNEX with default 0.5/0.5 weights when one basket genuinely different

- **Trigger:** Default mixture weights without elicitation.
- **Mechanism:** 50% EX weight still allows substantial borrowing.
- **Symptom:** Detected differential basket "softened" by borrowing.
- **Fix:** Sensitivity analysis over mixture weights (0.1, 0.3, 0.5, 0.7, 0.9); report range.

### Stratified randomisation ignored in subgroup analysis

- **Trigger:** Site-stratified randomisation; subgroup analysis without site adjustment.
- **Mechanism:** Achieved SE smaller than calculated SE.
- **Symptom:** Anti-conservative subgroup interaction tests.
- **Fix:** Include stratification factors as model covariates (Kahan-Morris 2012).

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| ~4x sample size for interaction detection vs main effect | Brookes et al 2004 *J Clin Epidemiol* 57:229 | Subgroup analyses are powered by interaction effect, not main effect |
| 0.5σ standardised effect = "noteworthy heterogeneity" | Dane et al 2019 EFSPI white paper *Pharm Stat* 18:126 | Discipline for declaring subgroup signals worth pursuing |
| Pre-specified subgroups for label claim | EMA 2019 subgroup guideline | Discovery subgroups cannot support indication restriction |
| Bonferroni: ~10 subgroups -> 30-50% power loss vs Hommel | Sarkar 1998 | Subgroup tests positively correlated; Bonferroni over-conservative |
| Honest splitting required for causal forest | Athey-Tibshirani-Wager 2019; Rehill 2025 | Without it, CIs are invalid; ~70% of papers omit |
| 11 credibility criteria | Sun BMJ 2012 344:e1553 | Canonical academic checklist; EMA 2019 implicitly references |
| RATE/AUTOC test (Yadlowsky 2025) | *JASA* 120(549):38-51 | Single-p omnibus for HTE-ranking predictive value |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| Per-subgroup p-values compared as evidence of HTE | Statistically invalid | Single model with interaction term; report interaction p, not per-subgroup p |
| "Non-significant interaction = no HTE" | Underpowered interaction test | Cite Brookes 2004 4x rule; "absence of evidence" framing |
| Breslow-Day non-sig taken as homogeneity | Low power with few strata | Forest plot + LR interaction; cite EMA 2019 |
| STEPP with naive pointwise CIs reported | Correlated estimates | Permutation supremum test (R `stepp::test_pattern`) |
| Causal forest CIs without honest splitting | Default settings | `honesty=TRUE`; tune; calibration test; cite Rehill 2025 |
| Bayesian shrinkage for "fishing" subgroup discovery | Misapplication | Use for replication planning per Hemmings-Koch; cite Dane white paper |
| Selected subgroup effect reported without shrinkage | Winner's curse | Apply Bayesian shrinkage (RBesT) or report selection-corrected estimate |
| Subgroup p < 0.05 claimed as evidence with no multiplicity | Common SAP omission | Graphical procedure with pre-specified alpha allocation (gMCP) |
| MH pooled OR with sign-reversing stratum ORs | Simpson's paradox | Report stratum-specific; switch to interaction model |
| Causal forest tau.hat reported without RATE test | Missing modern omnibus test | Yadlowsky 2025 RATE/AUTOC; cite |

## Anticipated Reviewer Pushback

| Pushback | Response |
|----------|----------|
| "Is this pre-specified?" | Per EMA 2019, distinguish assessment (pre-specified, regulatory weight) vs discovery (post-hoc, hypothesis-generating); be explicit. |
| "Where is multiplicity correction?" | Graphical procedure with alpha allocation per SAP; cite Bretz-Maurer 2009 + gMCP. |
| "Interaction power check?" | Brookes 2004: ~4x main effect SS needed; report interaction power calculation in SAP. |
| "Is this finding clinically plausible?" | Cite biological mechanism (e.g., known pharmacogenomic pathway); cite Sun 2012 criterion 11. |
| "Winner's curse control?" | Bayesian shrinkage per Dixon-Simon (RBesT) for adjusted estimate; or honest cross-validation. |
| "Why STEPP not categorical subgroups?" | Avoids arbitrary dichotomisation of continuous biomarker; cite Bonetti-Gelber 2004. |
| "Causal forest calibration?" | `test_calibration` p-value reported; honest splitting enabled; cite Rehill 2025 audit standards. |
| "Sensitivity to prior in shrinkage analysis?" | tau ~ HalfNormal(0, 0.1), 0.25, 0.5 reported; results stable across plausible priors. |

## References

- Athey S, Tibshirani J, Wager S. 2019. Generalized random forests. *Ann Stat* 47:1148-1178.
- Berry SM, Broglio KR, Groshen S, Berry DA. 2013. Bayesian hierarchical modeling of patient subpopulations: efficient designs of Phase II oncology clinical trials. *Clin Trials* 10:720-734.
- Bonetti M, Gelber RD. 2000. A graphical method to assess treatment-covariate interactions using the Cox model on subsets. *Biostatistics* 1:227-240.
- Bonetti M, Gelber RD. 2004. Patterns of treatment effects in subsets of patients in clinical trials. *Biostatistics* 5:465-481.
- Bretz F, Maurer W, Brannath W, Posch M. 2009. A graphical approach to sequentially rejective multiple test procedures. *Stat Med* 28:586-604.
- Brookes ST, Whitely E, Egger M, Davey Smith G, Mulheran PA, Peters TJ. 2004. Subgroup analyses in randomised controlled trials: quantifying the risks of false-positives and false-negatives. *J Clin Epidemiol* 57:229-236. (Author surname is "Whitely", not "Whitley".)
- Dane A, Spencer A, Rosenkranz G, Lipkovich I, Parke T. 2019. Subgroup analysis and interpretation for phase 3 confirmatory trials: white paper of the EFSPI/PSI working group. *Pharm Stat* 18:126-139.
- Dixon DO, Simon R. 1991. Bayesian subset analysis. *Biometrics* 47:871-881.
- Dusseldorp E, Van Mechelen I. 2014. Qualitative interaction trees: a tool to identify qualitative treatment-subgroup interactions. *Stat Med* 33:219-237.
- EMA. 2019. Guideline on the investigation of subgroups in confirmatory clinical trials. CHMP/EWP/117211/2010.
- Foster JC, Taylor JMG, Ruberg SJ. 2011. Subgroup identification from randomized clinical trial data. *Stat Med* 30:2867-2880.
- Gail M, Simon R. 1985. Testing for qualitative interactions between treatment effects and patient subsets. *Biometrics* 41:361-372.
- Goeman JJ, Hemerik J, Solari A. 2021. Only closed testing procedures are admissible for controlling false discovery proportions. *Ann Stat* 49:1218-1238.
- Hemmings R, Koch A. 2019. Commentary on Dane et al. *Pharm Stat* 18:140-141.
- Kahan BC, Morris TP. 2012. Improper analysis of trials randomised using stratified blocks or minimisation. *Stat Med* 31:328-340.
- Künzel SR, Sekhon JS, Bickel PJ, Yu B. 2019. Metalearners for estimating heterogeneous treatment effects using machine learning. *PNAS* 116:4156-4165.
- Lipkovich I, Dmitrienko A, Denne J, Enas G. 2011. Subgroup identification based on differential effect search: SIDES. *Stat Med* 30:2601-2621.
- Neuenschwander B, Wandel S, Roychoudhury S, Bailey S. 2016. Robust exchangeability designs for early phase clinical trials with multiple strata. *Pharm Stat* 15:123-134.
- Nie X, Wager S. 2021. Quasi-oracle estimation of heterogeneous treatment effects. *Biometrika* 108:299-319.
- Pan G, Wolfe DA. 1997. Test for qualitative interaction of clinical significance. *Biometrics* 53:1564-1570.
- Rehill A. 2025. A reliability review of applications of generalized random forests. *Int Stat Rev*.
- Schmidli H, Gsteiger S, Roychoudhury S, O'Hagan A, Spiegelhalter D, Neuenschwander B. 2014. Robust meta-analytic-predictive priors in clinical trials with historical control information. *Biometrics* 70:1023-1032.
- Senn S. 2018. Statistical pitfalls of personalised medicine. *Nature* 563:619-621.
- Sun X, Briel M, Walter SD, Guyatt GH. 2010. Is a subgroup effect believable? Updating criteria to evaluate the credibility of subgroup analyses. *BMJ* 340:c117.
- Sun X et al. 2012. Credibility of claims of subgroup effects in randomised controlled trials: systematic review. *BMJ* 344:e1553.
- Wager S, Athey S. 2018. Estimation and inference of heterogeneous treatment effects using random forests. *JASA* 113:1228-1242.
- Wang R, Lagakos SW, Ware JH, Hunter DJ, Drazen JM. 2007. Statistics in medicine -- reporting of subgroup analyses in clinical trials. *NEJM* 357:2189-2194.
- Weber S, Li Y, Seaman J, Kakizume T, Schmidli H. 2021. Applying meta-analytic-predictive priors with the R Bayesian evidence synthesis tools. *J Stat Softw* 100:19.
- Yadlowsky S, Fleming S, Shah N, Brunskill E, Wager S. 2025. Evaluating treatment prioritization rules via rank-weighted average treatment effects. *JASA* 120(549):38-51.

## Related Skills

- clinical-biostatistics/categorical-tests - CMH, Breslow-Day, chi-square within strata
- clinical-biostatistics/effect-measures - Forest plots, OR/RR/RD reporting
- clinical-biostatistics/logistic-regression - Interaction terms in regression
- clinical-biostatistics/multiplicity-graphical - Bretz-Maurer graphical procedures (in depth)
- clinical-biostatistics/bayesian-trials - MAP/EXNEX/Berry hierarchical models (in depth)
- clinical-biostatistics/trial-reporting - CONSORT 2025 + EMA 2019 reporting of subgroup analyses
- experimental-design/multiple-testing - General multiplicity correction methods
- machine-learning/biomarker-discovery - HTE for biomarker-defined subgroups

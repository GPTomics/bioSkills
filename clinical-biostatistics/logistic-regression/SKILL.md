---
name: bio-clinical-biostatistics-logistic-regression
description: Performs logistic regression for clinical trial outcomes including binary, ordinal, and multinomial models. Extracts odds ratios with confidence intervals, handles covariate adjustment, and provides Firth penalized regression for rare events or separation. Use when modeling binary or ordinal endpoints from clinical data.
tool_type: python
primary_tool: statsmodels
---

## Version Compatibility

Reference examples tested with: statsmodels 0.14+, scipy 1.12+, numpy 1.26+, pandas 2.1+, firthmodels 0.3+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Logistic Regression for Clinical Outcomes

**"Model clinical outcomes with logistic regression"** -> Fit generalized linear models with logit link to estimate treatment effects as odds ratios from clinical trial data.
- Python: `statsmodels.api.Logit()`, `statsmodels.formula.api.logit()`, `OrderedModel()`

## When to Use Logistic Regression vs Simpler Tests

Logistic regression answers a different question than a chi-square test: what is the adjusted treatment effect while controlling for covariates? Use chi-square or Fisher's exact when testing crude independence between treatment and outcome. Use logistic regression when confounders must be adjusted for, when continuous covariates (age, baseline severity) are part of the scientific model, or when the effect needs to be expressed as an OR. For RCTs with a simple two-arm design, balanced randomization, and no important prognostic covariates, a chi-square test is sufficient and more transparent.

## Binary Logistic Regression

**Goal:** Estimate the odds ratio for a treatment effect on a binary clinical outcome.

**Approach:** Fit a logistic regression model using statsmodels, extract coefficients, and exponentiate to obtain odds ratios with confidence intervals.

Three interfaces are available depending on the use case:

```python
import statsmodels.api as sm
import statsmodels.formula.api as smf
import numpy as np
import pandas as pd

# Formula API (recommended - auto-handles categoricals)
model = smf.logit('outcome ~ C(ARM, Treatment(reference="Placebo")) + age + C(sex)', data=df).fit()

# Matrix API
X = sm.add_constant(df[['treatment', 'age', 'sex_coded']])
model = sm.Logit(df['outcome'], X).fit()

# GLM (when a unified framework across link functions is needed)
model = sm.GLM(df['outcome'], X, family=sm.families.Binomial()).fit()
```

The formula API is preferred for clinical analyses because it automatically dummy-codes categorical variables and allows explicit reference level specification.

## Extracting Odds Ratios

**Goal:** Convert logistic regression coefficients to interpretable odds ratios with confidence intervals.

**Approach:** Exponentiate log-odds coefficients and their confidence bounds.

```python
or_table = pd.DataFrame({
    'OR': np.exp(model.params),
    'Lower_CI': np.exp(model.conf_int()[0]),
    'Upper_CI': np.exp(model.conf_int()[1]),
    'p_value': model.pvalues
})
or_table = or_table.drop('Intercept', errors='ignore')
```

An OR > 1 means the covariate increases the odds of the outcome. An OR < 1 means it decreases the odds. The confidence interval crossing 1.0 indicates non-significance at that alpha level.

### Clinical Significance vs Statistical Significance

A statistically significant OR does not imply clinical importance. With large trials (n > 5000), even trivially small effects reach p < 0.05. For many clinical outcomes, ORs between 0.8 and 1.25 are considered negligible regardless of p-value. Always interpret the OR magnitude and its confidence interval against pre-specified clinical thresholds (minimum clinically important difference), not just against the null value of 1.0. NNT provides a more clinically intuitive metric (see effect-measures skill).

## Interpreting Effect Direction

Correct interpretation depends entirely on how the reference category is set:

| Coding | Positive coefficient means |
|--------|--------------------------|
| treatment=1, placebo=0 | Treatment INCREASES odds |
| C(ARM, Treatment(reference='Placebo')) | Non-reference arm effect vs placebo |
| Alphabetical default C(ARM) | UNPREDICTABLE - depends on level order |

**CRITICAL:** statsmodels alphabetically orders categorical levels by default. 'Active' sorts before 'Placebo', which means Placebo is NOT the reference unless explicitly set. Always specify the reference: `C(ARM, Treatment(reference='Placebo'))`. Getting this wrong reverses the direction of the odds ratio.

```python
# WRONG: relies on alphabetical ordering
model = smf.logit('outcome ~ C(ARM)', data=df).fit()

# RIGHT: explicit reference level
model = smf.logit('outcome ~ C(ARM, Treatment(reference="Placebo"))', data=df).fit()
```

### Causal vs Associational Language

In RCTs, treatment effect ORs support causal language ("treatment reduced the odds") because randomization eliminates confounding. In observational studies, the same model yields only associational evidence ("treatment was associated with reduced odds"). Never use causal phrasing for observational data, even with covariate adjustment, because unmeasured confounding cannot be ruled out.

## Covariate Adjustment

**Goal:** Include pre-specified confounders to obtain unbiased treatment effect estimates.

**Approach:** Add confounders to the model formula based on clinical knowledge and pre-analysis plans, not data-driven selection.

| Variable type | Adjust? | Example |
|---------------|---------|---------|
| Confounder | Yes | Age, sex, baseline severity |
| Mediator | No - blocks causal path | Inflammation (if treatment affects it) |
| Collider | No - creates bias | Post-randomization variable caused by both treatment and outcome |

Use DAGs or domain knowledge to identify confounders. Stepwise selection inflates type I error and produces unstable models.

```python
# Pre-specified covariates from the statistical analysis plan
model = smf.logit(
    'outcome ~ C(ARM, Treatment(reference="Placebo")) + age + C(sex) + baseline_score',
    data=df
).fit()
```

### Effect Modification vs Confounding

A confounder distorts the treatment-outcome relationship and should be adjusted for. An effect modifier changes the magnitude of the treatment effect across its levels and should be tested via an interaction term, not adjusted away. Adjusting for an effect modifier without an interaction term forces a single average OR across subgroups, masking meaningful heterogeneity. If there is prior biological reason to believe a covariate modifies the treatment effect (e.g., a biomarker predicting drug response), include a treatment-by-covariate interaction term rather than just adding it as a main effect.

### Binary vs Ordinal Model Choice

When the outcome has ordered levels (e.g., none/mild/moderate/severe), the choice between dichotomizing to binary logistic and fitting an ordinal model is a genuine methodological decision:

| Factor | Favors binary | Favors ordinal |
|--------|--------------|----------------|
| Clinical question | "Does treatment prevent severe events?" (threshold question) | "Does treatment shift the entire severity distribution?" |
| Interpretability | Single OR for a clinically defined threshold, easy to communicate | Common OR across all cut-points, harder for clinical audiences |
| Statistical power | Lower -- discards ordering information | Higher -- preserves ordering (when proportional odds holds) |
| Proportional odds | Not required | Must hold; test before interpreting (see below) |
| Regulatory endpoint | Often pre-specified as binary (e.g., grade 3+ toxicity) | Less common in regulatory submissions |
| Category distribution | Most observations cluster in two groups | Events spread across multiple levels |

When the clinical question is inherently binary (e.g., severe toxicity yes/no), dichotomize even though ordinal is more "efficient." When the question is about the full distribution and the proportional odds assumption holds, ordinal preserves statistical power -- particularly valuable in small trials where discarding information is costly.

## Ordinal Logistic Regression

**Goal:** Model an ordinal outcome (e.g., mild/moderate/severe) while respecting category ordering.

**Approach:** Fit a proportional odds (cumulative logit) model using OrderedModel. Do not include an intercept -- thresholds serve this role.

```python
from statsmodels.miscmodels.ordinal_model import OrderedModel

df['severity'] = pd.Categorical(df['severity'], categories=['mild', 'moderate', 'severe'], ordered=True)

model = OrderedModel.from_formula('severity ~ treatment + age', data=df, distr='logit')
result = model.fit(method='bfgs')
print(result.summary())
```

OrderedModel internally removes the intercept because the threshold parameters (cut points between ordinal levels) replace it. Adding an explicit constant causes non-identifiability.

### Proportional Odds Assumption

The proportional odds model assumes the effect of each predictor is constant across all outcome cut-points. This must be tested -- if it fails, the model parameters are invalid. Compare the proportional odds model to a less constrained multinomial model via likelihood ratio test:

```python
from statsmodels.api import MNLogit

po_model = OrderedModel.from_formula('severity ~ treatment + age', data=df, distr='logit').fit(method='bfgs', disp=0)
mn_model = MNLogit.from_formula('severity ~ treatment + age', data=df).fit(disp=0)
lr_stat = 2 * (mn_model.llf - po_model.llf)
lr_df = mn_model.df_model - po_model.df_model
```

A significant LR test (p < 0.05) indicates the proportional odds assumption is violated. In that case, use multinomial logistic regression, or report stratum-specific ORs for each cut-point.

## Model Diagnostics

| Diagnostic | Method | Threshold |
|-----------|--------|-----------|
| Discrimination | ROC-AUC | > 0.7 acceptable, > 0.8 good |
| Pseudo R-squared | model.prsquared (McFadden) | > 0.2 excellent; not comparable to OLS R-squared |
| Calibration | Calibration plot (primary), Hosmer-Lemeshow (supplementary) | Predicted vs observed should follow diagonal |
| Prediction table | model.pred_table() | Rows=actual, columns=predicted |

### ROC-AUC

```python
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib.pyplot as plt

y_pred = model.predict()
auc = roc_auc_score(df['outcome'], y_pred)
fpr, tpr, _ = roc_curve(df['outcome'], y_pred)

fig, ax = plt.subplots()
ax.plot(fpr, tpr, label=f'AUC = {auc:.3f}')
ax.plot([0, 1], [0, 1], 'k--')
ax.set_xlabel('False positive rate')
ax.set_ylabel('True positive rate')
ax.legend()
plt.savefig('roc_curve.png', dpi=150)
```

### Hosmer-Lemeshow Test

**Goal:** Assess whether predicted probabilities match observed event rates across risk strata.

**Approach:** Bin subjects by predicted probability, compare observed vs expected counts per bin, and compute a chi-squared statistic.

```python
from scipy.stats import chi2

def hosmer_lemeshow(y_true, y_pred, n_groups=10):
    df_hl = pd.DataFrame({'y': y_true, 'prob': y_pred})
    df_hl['group'] = pd.qcut(df_hl['prob'], n_groups, duplicates='drop')
    grouped = df_hl.groupby('group').agg(obs=('y', 'sum'), n=('y', 'count'), pred=('prob', 'mean'))
    grouped['expected'] = grouped['n'] * grouped['pred']
    hl_stat = (((grouped['obs'] - grouped['expected']) ** 2) / (grouped['n'] * grouped['pred'] * (1 - grouped['pred']))).sum()
    actual_groups = len(grouped)
    p_value = 1 - chi2.cdf(hl_stat, actual_groups - 2)
    return hl_stat, p_value

hl_stat, hl_p = hosmer_lemeshow(df['outcome'], model.predict())
print(f'Hosmer-Lemeshow: chi2={hl_stat:.2f}, p={hl_p:.4f}')
```

The H-L test has low power in small samples (rarely rejects even for poor calibration when n < 200) and is oversensitive in large samples (rejects for trivial miscalibration when n > 2000). Use calibration plots (observed vs predicted rates in deciles) as the primary calibration assessment; H-L is supplementary. Note: `duplicates='drop'` can reduce the actual number of groups when predicted probabilities have ties, changing the degrees of freedom.

### Events Per Variable Rule

Logistic regression requires approximately 10 events per covariate for stable estimates. A model with 5 predictors (including treatment) needs at least 50 events. Below this threshold, coefficient estimates are biased away from zero, confidence intervals have poor coverage, and the model may not converge. Reduce the number of covariates or use Firth's method for sparse-event models.

### Modified Poisson for Common Outcomes

When the outcome is common (prevalence > 10%), OR substantially overestimates RR. Modified Poisson regression (log-Poisson with robust variance) directly estimates relative risks:

```python
import statsmodels.api as sm

poisson_model = sm.GLM(df['outcome'], sm.add_constant(df[['treatment', 'age']]),
                       family=sm.families.Poisson()).fit(cov_type='HC1')
rr = np.exp(poisson_model.params)
```

The `cov_type='HC1'` (Huber-White sandwich estimator) corrects the variance for the misspecified Poisson assumption. This approach is increasingly preferred in cohort studies and RCTs with common binary outcomes.

## Separation and Firth's Method

**Goal:** Handle complete or quasi-complete separation where standard MLE produces infinite coefficients.

**Approach:** Detect separation with a linear programming test, then apply Firth's penalized likelihood to obtain finite, bias-reduced estimates.

Detection signs: coefficients > 10, standard errors > 100, convergence warnings from the optimizer.

```python
from firthmodels import FirthLogisticRegression, detect_separation
import numpy as np

sep_result = detect_separation(X, y)
if sep_result.separation:
    print(sep_result.summary())

firth = FirthLogisticRegression()
firth.fit(X, y)
# firth.coef_ - coefficients (log-odds)
# firth.pvalues_ - Wald p-values
# firth.bse_ - standard errors
# firth.intercept_ - intercept term

or_firth = np.exp(firth.coef_)
```

Firth's method was originally designed for finite-sample bias reduction, not separation handling specifically. It adds a Jeffrey's prior penalty to the likelihood, which keeps coefficients finite even under separation. Also recommended for rare events (< 5% prevalence) where standard MLE has non-negligible bias.

Note: `firth.pvalues_` reports Wald p-values, which can be liberal with penalized estimates. The penalized likelihood ratio test (PLRT) is preferred for inference from Firth models. The firthmodels package provides PLRT via `firth.pvalues_lrt_` when available; otherwise compute manually by comparing penalized log-likelihoods of nested models.

For a statsmodels-compatible interface:

```python
from firthmodels.adapters.statsmodels import FirthLogit
import statsmodels.api as sm

X_const = sm.add_constant(X)
result = FirthLogit(y, X_const).fit()
print(result.summary())
# result.params, result.pvalues, result.conf_int()
```

## Common Pitfalls

- **Not setting reference category:** Alphabetical ordering silently reverses effect direction. Always use `C(ARM, Treatment(reference='Placebo'))`.
- **Using GLM when Logit suffices:** `sm.Logit` provides `pred_table()` and `prsquared` that `sm.GLM` does not. Use GLM only when switching between link functions.
- **Including intercept with OrderedModel:** Threshold parameters replace the intercept. An explicit constant causes non-identifiability and optimizer failures.
- **Adjusting for mediators:** Including post-treatment variables on the causal path between treatment and outcome attenuates the treatment effect toward the null.
- **Ignoring separation:** Standard MLE produces infinite coefficients with enormous standard errors. If any coefficient exceeds 10 or SE exceeds 100, check for separation and consider Firth's method.
- **Comparing pseudo R-squared to OLS:** McFadden's pseudo R-squared is not on the same scale as linear regression R-squared. Values above 0.2 indicate excellent fit.

## Related Skills

- clinical-biostatistics/cdisc-data-handling - Prepare analysis datasets from CDISC domains
- machine-learning/survival-analysis - Time-to-event modeling with Cox regression

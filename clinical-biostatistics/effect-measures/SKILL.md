---
name: bio-clinical-biostatistics-effect-measures
description: Computes and interprets treatment effect measures including odds ratios, risk ratios, number needed to treat, and confidence intervals from clinical trial data. Covers crude and adjusted measures, non-collapsibility of odds ratios, and forest plot visualization. Use when reporting treatment effects or comparing effect sizes across clinical studies.
tool_type: python
primary_tool: statsmodels
---

## Version Compatibility

Reference examples tested with: statsmodels 0.14+, numpy 1.26+, pandas 2.1+, matplotlib 3.8+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Treatment Effect Measures for Clinical Trials

**"Compute treatment effect sizes"** -> Calculate odds ratios, risk ratios, and number needed to treat with confidence intervals from clinical trial contingency tables or regression models.
- Python: `statsmodels.stats.contingency_tables.Table2x2()`, `np.exp(model.params)`

## Crude Measures from 2x2 Tables

**Goal:** Compute unadjusted odds ratios and risk ratios directly from a contingency table.

**Approach:** Use Table2x2 for comprehensive 2x2 analysis including OR, RR, and their confidence intervals.

```python
from statsmodels.stats.contingency_tables import Table2x2
import numpy as np

# Table layout: [[exposed+outcome, exposed+no_outcome], [unexposed+outcome, unexposed+no_outcome]]
table = np.array([[a, b], [c, d]])
t = Table2x2(table)
print(t.oddsratio)
print(t.oddsratio_confint())
print(t.log_oddsratio_se)
print(t.riskratio)
print(t.riskratio_confint())
```

**Table orientation is critical**: Table2x2 interprets the first column as "outcome present." When using `pd.crosstab()`, ensure the event column (1) is ordered first: `cross = cross[[1, 0]]`. If columns are reversed, the OR is the reciprocal of the intended value -- a silent, direction-reversing error.

`oddsratio_confint()` and `riskratio_confint()` return Wald-type confidence intervals by default. For the log scale, `log_oddsratio_confint()` and `log_riskratio_confint()` are also available. All accept an `alpha` parameter (default 0.05 for 95% CI).

## Choosing the Right Measure

| Measure | Study design | Outcome | When to use |
|---------|-------------|---------|-------------|
| Odds Ratio | Case-control, cross-sectional, RCT | Binary (cumulative) | Default for logistic regression; only valid measure in case-control |
| Risk Ratio | Cohort, RCT | Binary (cumulative) | Preferred for common outcomes (>10% prevalence) in cohort/RCT |
| Hazard Ratio | Any with time-to-event | Survival | When follow-up varies or time-to-event matters |

When outcome prevalence exceeds 10%, the OR overestimates the RR. An OR of 2.0 with 30% baseline risk corresponds to RR of approximately 1.5. ORs should not be interpreted as if they were RRs in this regime.

### Reporting ORs as Percentage Changes

Clinical reports often express ORs as percentage changes in odds: `(OR - 1) * 100%`. An OR of 1.524 means "increases odds by 52.4%"; an OR of 0.755 means "decreases odds by 24.5%." Note the asymmetry: OR 0.5 is a 50% reduction, but its reciprocal OR 2.0 is a 100% increase.

This percentage refers to the change in *odds*, not the change in *risk* (probability). Clinicians and lay audiences routinely conflate these. When baseline risk is low (< 10%), OR approximates RR and the percentage translation roughly holds for risk too. For common outcomes, the percentage change in odds substantially overstates the percentage change in risk. Always clarify "percentage change in odds" when reporting, or convert to risk difference or NNT (see below) for clinical audiences.

### Clinical Significance Benchmarks

Effect measures must be interpreted against clinical context, not just statistical thresholds. General benchmarks: OR 1.0-1.5 (small effect, may lack clinical relevance unless the outcome is severe), OR 1.5-3.0 (moderate, often clinically meaningful), OR > 3.0 (large). However, these are crude -- an OR of 1.2 for mortality may be more important than an OR of 3.0 for mild headache. Always evaluate the OR against the severity and reversibility of the outcome, the baseline risk, and any pre-specified minimum clinically important difference (MCID).

### Causal Interpretation by Study Design

The same OR supports different levels of causal inference depending on study design. In an RCT, OR = 0.6 for treatment vs placebo can be stated as "treatment reduced the odds." In a cohort study with adequate adjustment, the same OR supports "treatment was associated with reduced odds." In a case-control study, only the OR is estimable (not RR), and it reflects association under the assumption that controls represent the source population. Always pair the effect estimate with the study design when reporting.

## Adjusted OR from Logistic Regression

**Goal:** Extract covariate-adjusted odds ratios from a fitted logistic regression model.

**Approach:** Exponentiate the model coefficients and confidence interval bounds.

```python
import statsmodels.formula.api as smf
import pandas as pd

model = smf.logit('outcome ~ C(ARM, Treatment(reference="Placebo")) + age + C(sex)', data=df).fit()
or_table = pd.DataFrame({
    'OR': np.exp(model.params),
    'Lower_CI': np.exp(model.conf_int()[0]),
    'Upper_CI': np.exp(model.conf_int()[1]),
    'p_value': model.pvalues
})
```

The exponentiated intercept is not clinically meaningful and is typically dropped from the results table.

## Confidence Interval Methods

| Method | When to use | Properties |
|--------|-------------|-----------|
| Wald (default) | Large samples (n > 100) | Fast, symmetric on log scale; poor coverage for small n |
| Profile likelihood | Small samples (n < 50), near boundaries | Transformation-invariant, asymmetric, more accurate |
| Bootstrap | Complex models, clustered data | No distributional assumptions; computationally expensive |

Wald CIs are the default in virtually all statistical software but have the poorest coverage properties for small samples. Profile likelihood CIs are superior because they are transformation-invariant -- the CI for OR equals the exponentiated CI for log(OR). For critical analyses with small samples, profile likelihood is recommended.

## Number Needed to Treat

**Goal:** Convert effect measures to a clinically intuitive metric -- the number of patients who must be treated to prevent one additional adverse event (NNT) or cause one additional harm (NNH).

**Approach:** Compute the absolute risk reduction from event rates, then invert.

```python
def nnt_from_rates(treated_events, treated_n, control_events, control_n):
    eer = treated_events / treated_n
    cer = control_events / control_n
    arr = abs(cer - eer)
    return int(np.ceil(1 / arr)) if arr > 0 else float('inf')

def nnt_from_or(odds_ratio, baseline_risk):
    baseline_odds = baseline_risk / (1 - baseline_risk)
    treatment_odds = baseline_odds * odds_ratio
    treatment_risk = treatment_odds / (1 + treatment_odds)
    arr = abs(baseline_risk - treatment_risk)
    return int(np.ceil(1 / arr)) if arr > 0 else float('inf')
```

NNT computed from an OR requires knowing the baseline risk. The same OR yields very different NNTs at different baseline rates:

| OR | Baseline risk 5% | Baseline risk 20% | Baseline risk 50% |
|----|-------------------|--------------------|--------------------|
| 0.5 | NNT = 42 | NNT = 12 | NNT = 6 |
| 0.7 | NNT = 70 | NNT = 20 | NNT = 12 |

Always report baseline risk alongside NNT. An NNT without context is uninterpretable. NNT confidence intervals are non-trivial: when the ARR CI crosses zero, the NNT CI passes through infinity (NNTB -> infinity -> NNTH). Use the Altman-Andersen method or bootstrap for NNT CIs.

### Risk Difference

The risk difference (RD = p_treated - p_control) is the most directly interpretable absolute measure and is increasingly emphasized by regulators alongside OR/RR. Table2x2 does not provide RD directly; compute manually:

```python
p_treated = a / (a + b)
p_control = c / (c + d)
rd = p_treated - p_control
se_rd = np.sqrt(p_treated * (1 - p_treated) / (a + b) + p_control * (1 - p_control) / (c + d))
rd_ci = (rd - 1.96 * se_rd, rd + 1.96 * se_rd)
```

RD is collapsible (unlike OR): the marginal and conditional RDs are equal absent confounding, making it straightforward to compare across adjusted and unadjusted analyses.

## Non-Collapsibility of the Odds Ratio

A conditional (adjusted) OR typically differs from the marginal (unadjusted) OR even without confounding. Conditioning on a prognostic covariate -- one that predicts outcome but is unrelated to treatment -- changes the OR. In most clinical settings, the adjusted OR tends to be further from the null than the marginal OR, though the direction depends on the covariate-outcome correlation structure. This is a mathematical property of the OR, not bias.

Risk ratios and risk differences do not have this property. They are collapsible: the marginal and conditional measures are equal in the absence of confounding.

Practical implication: when comparing ORs across studies or between adjusted and unadjusted models, differences may reflect non-collapsibility rather than confounding. This distinction matters for meta-analysis and regulatory submissions where adjusted and unadjusted ORs are both reported.

## Forest Plots

**Goal:** Visualize treatment effects across subgroups or studies on a common scale.

**Approach:** Plot point estimates with confidence interval error bars on a log-scaled x-axis with a reference line at OR = 1.

```python
import matplotlib.pyplot as plt

def forest_plot(labels, odds_ratios, lower_cis, upper_cis, figsize=(8, 6)):
    fig, ax = plt.subplots(figsize=figsize)
    y_pos = range(len(labels))
    ax.errorbar(odds_ratios, y_pos,
                xerr=[np.array(odds_ratios) - np.array(lower_cis),
                      np.array(upper_cis) - np.array(odds_ratios)],
                fmt='D', color='black', capsize=3, markersize=5)
    ax.axvline(x=1.0, color='gray', linestyle='--', linewidth=0.8)
    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(labels)
    ax.set_xlabel('Odds Ratio (95% CI)')
    ax.set_xscale('log')
    plt.tight_layout()
    return fig
```

The log scale ensures symmetric visual representation of reciprocal effects (OR 0.5 and OR 2.0 are equidistant from OR 1.0). The dashed vertical line at 1.0 represents no effect; CIs crossing this line indicate non-significance.

## Common Pitfalls

- Interpreting OR as RR when prevalence > 10%: systematically overstates the treatment effect
- Comparing adjusted and unadjusted ORs to quantify confounding: non-collapsibility contributes to the difference
- Computing NNT from OR without specifying baseline risk: the same OR gives very different NNTs
- Using Wald CIs with small samples: poor coverage; prefer profile likelihood
- Reporting OR without CI: effect size without uncertainty is uninterpretable
- Plotting forest plots on a linear scale: reciprocal effects appear asymmetric, distorting visual interpretation

## Related Skills

- clinical-biostatistics/logistic-regression - Fit models that produce adjusted ORs
- clinical-biostatistics/categorical-tests - Contingency table tests that complement effect measures
- clinical-biostatistics/subgroup-analysis - Forest plots for subgroup effects
- clinical-biostatistics/trial-reporting - CONSORT-compliant effect reporting
- machine-learning/survival-analysis - Hazard ratios for time-to-event outcomes

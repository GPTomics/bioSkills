---
name: bio-clinical-biostatistics-subgroup-analysis
description: Performs stratified and subgroup analyses for clinical trial data. Covers Mantel-Haenszel pooling, Breslow-Day homogeneity testing, interaction terms in regression, multiple comparisons correction, and forest plot visualization. Use when analyzing treatment effects across patient subgroups or controlling for stratification variables.
tool_type: python
primary_tool: statsmodels
---

## Version Compatibility

Reference examples tested with: statsmodels 0.14+, scipy 1.12+, numpy 1.26+, pandas 2.1+, matplotlib 3.8+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Subgroup Analysis

**"Analyze treatment effects across subgroups"** -> Test whether treatment effects differ across patient subgroups using stratified analysis, interaction terms, and multiplicity-adjusted comparisons.
- Python: `statsmodels.stats.contingency_tables.StratifiedTable()`, `statsmodels.formula.api.logit()`

## Mantel-Haenszel Stratified Analysis

**Goal:** Estimate a pooled treatment effect across strata while allowing for different baseline rates.

**Approach:** Construct per-stratum 2x2 tables and compute the Mantel-Haenszel weighted odds ratio.

```python
from statsmodels.stats.contingency_tables import StratifiedTable
import pandas as pd
import numpy as np

tables = []
for stratum in df['subgroup'].unique():
    sub = df[df['subgroup'] == stratum]
    t = pd.crosstab(sub['treatment'], sub['outcome']).values
    if t.shape == (2, 2):
        tables.append(t)

st = StratifiedTable(tables)
print(st.summary())
print(st.oddsratio_pooled)              # MH pooled OR
print(st.oddsratio_pooled_confint())    # 95% CI
print(st.test_null_odds())              # H0: common OR = 1
print(st.test_equal_odds())             # Breslow-Day: H0: all stratum ORs equal
```

## Breslow-Day Test for Homogeneity

The `test_equal_odds()` method tests whether stratum-specific odds ratios are equal. A significant result (p < 0.05) suggests effect modification across strata.

Breslow-Day has low power with few strata or small stratum sizes. A non-significant result does not prove homogeneity -- it may reflect insufficient power to detect heterogeneity. With many strata, visual assessment via forest plot should supplement the formal test before drawing conclusions.

## Interaction Terms in Regression

**Goal:** Test whether the treatment effect varies by a subgroup variable.

**Approach:** Fit a single logistic model with an interaction term, where the interaction coefficient tests effect modification directly.

```python
import statsmodels.formula.api as smf

# Single model with interaction -- the correct approach
model = smf.logit(
    'outcome ~ C(treatment, Treatment(reference="Placebo")) * C(age_group)', data=df
).fit()
# The interaction coefficient tests whether treatment effect differs by age group

# Extract subgroup-specific ORs for reporting
for group in df['age_group'].unique():
    sub_model = smf.logit(
        'outcome ~ C(treatment, Treatment(reference="Placebo"))',
        data=df[df['age_group'] == group]
    ).fit()
    or_val = np.exp(sub_model.params.iloc[1])
    ci = np.exp(sub_model.conf_int().iloc[1])
    print(f'{group}: OR={or_val:.3f} ({ci[0]:.3f}-{ci[1]:.3f})')
```

The proper way to test for subgroup effects is via interaction terms in a single model, NOT by comparing p-values from separate per-subgroup models. Separate models have different power and comparing their p-values is statistically invalid.

### When to Suspect Effect Modification

Not every baseline variable warrants an interaction test. Variables should be tested for effect modification when there is prior scientific reason to expect the treatment effect to differ across levels. Common biologically motivated effect modifiers include: disease severity (treatment may work only in severe disease), genetic variants affecting drug metabolism (pharmacogenomic subgroups), biomarkers of the targeted pathway, and age when pharmacokinetics differ substantially. Testing every available demographic variable without scientific rationale inflates false positives and produces uninterpretable results, even with multiplicity correction.

## Multiplicative vs Additive Interaction

Logistic regression tests multiplicative interaction (ratio of ORs). Null multiplicative interaction does not imply null additive interaction. For public health decisions, additive interaction is often more relevant.

RERI (Relative Excess Risk due to Interaction) measures additive interaction on the multiplicative scale:

```python
# From a model: outcome ~ treatment + subgroup_indicator + treatment:subgroup_indicator
# OR_11 = OR for treated subjects in the subgroup
# OR_10 = OR for treated subjects not in the subgroup
# OR_01 = OR for untreated subjects in the subgroup
reri = or_11 - or_10 - or_01 + 1
```

RERI = 0 indicates no additive interaction. Positive RERI indicates synergism (combined effect exceeds sum of individual effects). RERI CIs require the delta method or bootstrap, since RERI is a nonlinear function of the ORs.

### Quantitative vs Qualitative Interaction

Quantitative interaction means the treatment effect varies in magnitude across subgroups but remains in the same direction. Qualitative (crossover) interaction means the effect reverses direction -- beneficial in one subgroup, harmful in another. This distinction is critical for regulatory decisions: qualitative interaction may warrant restricting the indication to the benefiting subgroup. The Gail-Simon test formally tests for qualitative interaction.

### Power for Interaction Detection

Detecting a treatment-by-subgroup interaction requires approximately 4 times the sample size needed to detect the main treatment effect. A trial powered to detect OR = 0.6 overall cannot reliably detect subgroup differences of similar magnitude. Non-significant interaction tests should be interpreted cautiously -- absence of evidence is not evidence of absence, especially in underpowered subgroup analyses.

## Multiple Comparisons in Subgroup Analyses

**Goal:** Control error rates when testing treatment effects across multiple subgroups.

**Approach:** Apply FWER or FDR correction to the set of subgroup-specific p-values.

```python
from statsmodels.stats.multitest import multipletests

subgroup_pvalues = [0.03, 0.15, 0.04, 0.22, 0.01]

# FWER control (appropriate for regulatory/confirmatory)
reject_fwer, adjusted_fwer, _, _ = multipletests(subgroup_pvalues, method='holm')

# FDR control (appropriate for exploratory)
reject_fdr, adjusted_fdr, _, _ = multipletests(subgroup_pvalues, method='fdr_bh')
```

| Method | Controls | Use case |
|--------|---------|----------|
| Holm (step-down Bonferroni) | FWER | Confirmatory/regulatory subgroup tests |
| Hochberg (step-up) | FWER | Less conservative than Holm; valid only under independence or positive regression dependency (PRDS) |
| Benjamini-Hochberg | FDR | Exploratory subgroup screening |

## Pre-Specified vs Post-Hoc Subgroups

| Aspect | Pre-specified | Post-hoc |
|--------|--------------|----------|
| Timing | Before unblinding, in SAP | After seeing data |
| Credibility | High (if biologically justified) | Low (hypothesis-generating only) |
| Regulatory weight | Can support labeling claims | Cannot support claims alone |
| Multiplicity adjustment | Required per SAP | Required + heavy skepticism |

EMA credibility criteria for subgroup claims: (1) pre-specified and biologically plausible, (2) significant interaction test, (3) consistent across related endpoints, (4) ideally replicated in an independent study.

### Evaluating Biological Plausibility

A statistically significant interaction does not establish a real subgroup effect. Before concluding that treatment works differently in a subgroup, ask: (1) Is there a known biological mechanism? A drug targeting estrogen receptors showing differential efficacy by sex has mechanistic support. An age subgroup effect for an antibiotic does not. (2) Is the direction consistent with the mechanism? (3) Are related biomarkers concordant? Absent biological rationale, a "significant" subgroup finding among many tested subgroups is more likely a false positive than a real effect, regardless of the p-value.

## Forest Plots for Subgroup Effects

**Goal:** Visualize point estimates and confidence intervals across subgroups on a common scale.

**Approach:** Plot subgroup-specific ORs with error bars on a log-scaled axis with a reference line at OR = 1.

```python
import matplotlib.pyplot as plt
import numpy as np

def subgroup_forest_plot(labels, ors, lower_cis, upper_cis, overall_or=None, figsize=(8, 6)):
    fig, ax = plt.subplots(figsize=figsize)
    y_pos = range(len(labels))
    ax.errorbar(ors, y_pos,
                xerr=[np.array(ors) - np.array(lower_cis),
                      np.array(upper_cis) - np.array(ors)],
                fmt='D', color='black', capsize=3, markersize=5)
    ax.axvline(x=1.0, color='gray', linestyle='--', linewidth=0.8)
    if overall_or is not None:
        ax.axvline(x=overall_or, color='blue', linestyle=':', linewidth=0.8, alpha=0.5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels)
    ax.set_xlabel('Odds Ratio (95% CI)')
    ax.set_xscale('log')
    plt.tight_layout()
    return fig
```

## Common Pitfalls

- **Comparing per-subgroup p-values:** Invalid; separate models have different power and comparing their p-values conflates effect size with sample size. Always use an interaction term in a single model.
- **Non-significant = no difference:** Absence of evidence is not evidence of absence, especially with low power in small subgroups.
- **Post-hoc fishing as confirmatory:** Regulatory bodies reject subgroup claims that were not pre-specified in the statistical analysis plan.
- **Breslow-Day power:** Non-significance with few strata does not confirm homogeneity. Supplement with visual assessment.
- **Multiplicity ignorance:** Testing 10 subgroups at alpha=0.05 yields ~40% probability of at least one false positive under the global null.

## Related Skills

- clinical-biostatistics/categorical-tests - Chi-square and CMH tests used within strata
- clinical-biostatistics/effect-measures - OR computation and forest plots
- clinical-biostatistics/logistic-regression - Interaction terms in regression models
- clinical-biostatistics/trial-reporting - CONSORT-compliant subgroup reporting
- experimental-design/multiple-testing - General multiplicity correction methods

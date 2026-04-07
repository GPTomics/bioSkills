---
name: bio-clinical-biostatistics-categorical-tests
description: Tests associations between categorical variables in clinical data using chi-square, Fisher's exact, and Cochran-Mantel-Haenszel tests. Computes effect sizes and post-hoc pairwise comparisons. Use when analyzing categorical outcomes or testing treatment-outcome independence in clinical trials.
tool_type: python
primary_tool: scipy
---

## Version Compatibility

Reference examples tested with: scipy 1.12+, statsmodels 0.14+, pingouin 0.5+, pandas 2.1+, numpy 1.26+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Categorical Association Tests for Clinical Data

**"Test association between categorical variables"** -> Determine whether treatment assignment and a categorical clinical outcome are statistically independent using contingency table tests.
- Python: `scipy.stats.chi2_contingency()`, `scipy.stats.fisher_exact()`, `statsmodels.stats.contingency_tables.StratifiedTable()`

## Test Selection

| Criterion | Chi-square (Pearson) | Fisher's exact |
|-----------|---------------------|----------------|
| Expected cell counts | All >= 5 | Any (especially < 5) |
| Table size | Any RxC | 2x2 only via `scipy.stats.fisher_exact()`; for RxC with low expected counts, use chi-square or a permutation test |
| Sample size | Large (n > 40 for 2x2) | Any |
| P-value type | Asymptotic | Exact |

Decision flowchart:

```
Binary outcome, two groups:
  Expected cells >= 5 --> Pearson chi-square (correction=False)
  Expected cells < 5  --> Fisher's exact
  Stratified by third variable --> Cochran-Mantel-Haenszel
  Paired/matched data --> McNemar's test
```

### Matching Test to Scientific Question

Beyond sample size criteria, the choice of test reflects the scientific question. A chi-square test asks "are treatment and outcome independent?" -- a symmetric question with no directionality or magnitude. When the goal is to quantify how much treatment changes the odds of outcome, Fisher's exact (which returns an OR) or logistic regression is more appropriate. CMH is not merely "stratified chi-square" -- it asks whether a common treatment effect exists across strata, which is a question about generalizability. McNemar's addresses a longitudinal question: did subjects' status change? Choose the test that answers the question being asked, not just the one whose assumptions are satisfied.

## Chi-Square Test

**Goal:** Test whether treatment group and outcome category are independent.

**Approach:** Construct a contingency table and compute the Pearson chi-square statistic.

```python
from scipy.stats import chi2_contingency, fisher_exact
import pandas as pd

table = pd.crosstab(df['treatment'], df['outcome'])
chi2, p, dof, expected = chi2_contingency(table, correction=False)
```

The `correction=True` default in scipy applies Yates' continuity correction for 2x2 tables. Current consensus: Yates' correction is overly conservative and inflates Type II error. Standard practice is `correction=False` for Pearson chi-square. For small samples where asymptotic approximation breaks down, Fisher's exact test is the appropriate alternative rather than a continuity correction.

## Fisher's Exact Test

**Goal:** Test association when expected cell counts are too small for chi-square approximation.

**Approach:** Compute an exact p-value from the hypergeometric distribution.

```python
odds_ratio, p_fisher = fisher_exact(table.values, alternative='two-sided')
```

Since scipy 1.10, the returned `odds_ratio` is the sample (unconditional) odds ratio, not the conditional MLE. For the conditional MLE matching R's `fisher.test`, use `scipy.stats.contingency.odds_ratio(table, kind='conditional')`. Fisher's exact test is appropriate whenever any expected cell count falls below 5, or when the total sample size is small (n < 40 for 2x2 tables).

## Effect Sizes

**Goal:** Quantify the strength of association beyond statistical significance.

**Approach:** Compute phi coefficient (2x2) or Cramer's V (RxC) from the chi-square statistic.

```python
import numpy as np

n = table.values.sum()
phi = np.sqrt(chi2 / n)

k = min(table.shape) - 1
cramers_v = np.sqrt(chi2 / (n * k))
```

Effect size benchmarks (Cohen, 1988):

| df | Small | Medium | Large |
|----|-------|--------|-------|
| 1 | 0.10 | 0.30 | 0.50 |
| 2 | 0.07 | 0.21 | 0.35 |
| 3 | 0.06 | 0.17 | 0.29 |

Phi is equivalent to Cramer's V for 2x2 tables (where k = 1). For larger tables, only Cramer's V is valid because phi can exceed 1.0.

Bias-corrected Cramer's V via pingouin (corrects upward bias in small samples):

```python
import pingouin as pg

expected, observed, stats = pg.chi2_independence(df, x='treatment', y='outcome')
```

The `stats` DataFrame contains columns for `test`, `lambda`, `chi2`, `dof`, `pval`, `cramer`, and `power` across multiple test variants (Pearson, Cressie-Read, log-likelihood, Freeman-Tukey, mod-log-likelihood, Neyman).

## Post-Hoc Pairwise Comparisons

**Goal:** Identify which specific category pairs drive a significant overall chi-square result.

**Approach:** Perform chi-square tests on all pairwise subsets and adjust p-values for multiple comparisons.

```python
from statsmodels.stats.multitest import multipletests
from itertools import combinations

categories = df['outcome'].unique()
pvalues = []
comparisons = []
for cat1, cat2 in combinations(categories, 2):
    subset = df[df['outcome'].isin([cat1, cat2])]
    sub_table = pd.crosstab(subset['treatment'], subset['outcome'])
    _, p_val, _, _ = chi2_contingency(sub_table, correction=False)
    pvalues.append(p_val)
    comparisons.append(f'{cat1} vs {cat2}')

reject, adjusted_p, _, _ = multipletests(pvalues, method='holm')
```

The `multipletests` `method` parameter controls the correction procedure. Common choices:
- `'bonferroni'`: most conservative, single-step FWER control
- `'holm'`: step-down Bonferroni, less conservative FWER control (note: the actual default is `'hs'` / Holm-Sidak; always specify the method explicitly)
- `'fdr_bh'`: Benjamini-Hochberg FDR, appropriate for exploratory analyses

For regulatory submissions, FWER methods (`'bonferroni'`, `'holm'`) are standard. For exploratory biomarker analyses, FDR (`'fdr_bh'`) may be acceptable.

## Cochran-Mantel-Haenszel Test

**Goal:** Test treatment-outcome association while controlling for a stratification variable.

**Approach:** Construct a series of 2x2 tables (one per stratum) and compute the Mantel-Haenszel pooled odds ratio and test statistic.

```python
from statsmodels.stats.contingency_tables import StratifiedTable

tables = []
for stratum in df['site'].unique():
    stratum_data = df[df['site'] == stratum]
    t = pd.crosstab(stratum_data['treatment'], stratum_data['outcome']).values
    if t.shape == (2, 2):
        tables.append(t)

st = StratifiedTable(tables)
print(st.test_null_odds())
print(st.oddsratio_pooled)
print(st.oddsratio_pooled_confint())
```

`test_null_odds()` is the Mantel-Haenszel test of H0: common OR = 1 across all strata. `oddsratio_pooled` is the MH pooled estimate. `test_equal_odds()` is the Breslow-Day test for homogeneity of odds ratios across strata -- a significant result indicates the OR varies by stratum, and the pooled OR is not a valid summary; report stratum-specific ORs instead.

CMH assumes no three-way interaction between treatment, outcome, and the stratifying variable. It also assumes the treatment-outcome association is the same direction across strata. If the association reverses direction across strata (Simpson's paradox), the pooled MH OR can be misleading even when homogeneity holds on average. Always examine stratum-specific ORs alongside the pooled estimate to detect sign reversals.

## McNemar's Test

For paired categorical data (e.g., pre/post treatment on the same subjects):

```python
from statsmodels.stats.contingency_tables import mcnemar

# table[i,j] = count with outcome i at time 1 and j at time 2
table = np.array([[45, 15], [5, 35]])
result = mcnemar(table, exact=True)
print(result.statistic, result.pvalue)
```

The `exact=True` parameter uses the binomial distribution for exact p-values; appropriate when the number of discordant pairs is small (< 25). With `exact=False`, a chi-square approximation is used, and the `correction` parameter controls Yates' continuity correction.

## Common Pitfalls

- Using `correction=True` (Yates') by default: overly conservative; use Fisher's exact for small samples instead
- Comparing chi-square p-values across subgroups instead of using CMH or interaction tests
- Ignoring expected cell counts: chi-square approximation is unreliable when expected counts < 5
- Using chi-square on tables with structural zero cells: switch to Fisher's exact (2x2) or exact permutation test; do not add continuity corrections to chi-square for zero cells
- Applying FDR correction to confirmatory analyses: FWER methods (Bonferroni, Holm) are appropriate for regulatory submissions
- Computing phi on tables larger than 2x2: phi can exceed 1.0; use Cramer's V for RxC tables
- Ignoring study design in interpretation: a significant chi-square in an RCT supports a causal treatment effect because randomization controls confounding; the same result in an observational study indicates association only and needs CMH or logistic regression with covariate adjustment to approach causal inference

## Related Skills

- clinical-biostatistics/effect-measures - Detailed OR and risk ratio analysis
- clinical-biostatistics/logistic-regression - Regression-based alternative to contingency tests
- clinical-biostatistics/subgroup-analysis - Stratified analysis with interaction terms
- experimental-design/multiple-testing - General multiple testing correction methods

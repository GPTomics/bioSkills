---
name: bio-clinical-biostatistics-trial-reporting
description: Prepares statistical tables and reports for clinical trials following regulatory standards. Generates Table 1 baseline characteristics, defines analysis populations (ITT, per-protocol, safety), performs multiple imputation for missing data, and follows CONSORT and ICH E9 guidelines. Use when creating analysis reports, handling missing data, or preparing regulatory submissions from clinical trials.
tool_type: python
primary_tool: tableone
---

## Version Compatibility

Reference examples tested with: tableone 0.9+, statsmodels 0.14+, sklearn 1.4+, pandas 2.1+, numpy 1.26+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Trial Reporting

**"Prepare clinical trial statistical reports"** -> Generate baseline characteristics tables, define analysis populations, handle missing data through multiple imputation, and structure results per CONSORT/ICH E9 standards.
- Python: `tableone.TableOne()`, `sklearn.impute.IterativeImputer()`

## Table 1 Baseline Characteristics

**Goal:** Summarize and compare baseline demographics and clinical characteristics across treatment arms.

**Approach:** Use TableOne to generate a publication-ready baseline table with optional p-values and standardized mean differences.

```python
from tableone import TableOne

columns = ['age', 'sex', 'race', 'bmi', 'baseline_score', 'disease_stage']
categorical = ['sex', 'race', 'disease_stage']

table1 = TableOne(df, columns=columns, categorical=categorical,
                  groupby='ARM', pval=True, smd=True,
                  missing=True, overall=True)
print(table1.tabulate(tablefmt='github'))
table1.to_excel('table1.xlsx')
```

In RCTs, baseline imbalance is due to chance by definition. Significance tests in Table 1 test whether randomization produced balanced groups, but randomization is a known mechanism, not a hypothesis. Standardized mean differences (SMD > 0.1 suggests meaningful imbalance) are more informative than p-values for assessing balance. CONSORT 2010 discouraged baseline p-values; many journals still require them.

## Analysis Populations

| Population | Definition | Bias direction | Primary use |
|-----------|-----------|---------------|-------------|
| ITT (Full Analysis Set) | All randomized, as randomized | Conservative (toward null) | Primary analysis (regulatory standard) |
| Per-Protocol | Completed treatment per protocol | Anti-conservative (inflates effect) | Sensitivity analysis |
| Modified ITT | ITT excluding never-treated | Middle ground | Common in practice |
| Safety | All received at least one dose | N/A | Adverse event analysis |

ITT preserves the benefits of randomization and prevents bias from differential dropout. ICH E9 recommends ITT as the primary analysis population.

```python
# ITT: all randomized subjects
itt = dm.copy()

# Per-protocol: completed treatment without major violations
pp = dm[dm['USUBJID'].isin(completers) & ~dm['USUBJID'].isin(protocol_violators)]

# Safety: received at least one dose
dosed = ex[ex['EXDOSE'] > 0]['USUBJID'].unique()
safety = dm[dm['USUBJID'].isin(dosed)]

for name, pop in [('ITT', itt), ('Per-Protocol', pp), ('Safety', safety)]:
    print(f'{name}: n={len(pop)}, arms={pop["ARM"].value_counts().to_dict()}')
```

## Missing Data Mechanisms

| Mechanism | Definition | Testable? | Valid method |
|-----------|-----------|-----------|-------------|
| MCAR | Independent of all data | Partially (Little's test) | Complete-case unbiased but loses power |
| MAR | Depends on observed data only | No (assumption) | Multiple imputation valid |
| MNAR | Depends on unobserved values | No | Requires sensitivity analysis |

MAR vs MNAR cannot be distinguished from observed data alone. This is a fundamental limitation. The assumed mechanism should be pre-specified in the statistical analysis plan.

### Reasoning About Why Data Is Missing

The MCAR/MAR/MNAR classification is a statistical abstraction. Clinical reasoning requires examining the actual reasons patients have missing data. Did the patient drop out due to adverse events (likely MNAR -- the missing outcome is related to the outcome itself)? Did a site close early for administrative reasons (likely MCAR)? Did sicker patients miss follow-up visits (MAR if sickness is captured in observed covariates, MNAR if not)? Examining the DS (Disposition) domain to tabulate reasons for discontinuation by treatment arm is essential. If discontinuation rates or reasons differ between arms, missing data is likely informative and sensitivity analyses under MNAR are mandatory.

## Multiple Imputation with Rubin's Rules

**Goal:** Generate multiple plausible complete datasets and pool results to properly account for imputation uncertainty.

**Approach:** Use IterativeImputer with `sample_posterior=True` to generate m imputed datasets, fit the analysis model on each, and combine estimates using Rubin's rules.

```python
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
import statsmodels.formula.api as smf
import numpy as np
import pandas as pd

n_imputations = 20
imputer = IterativeImputer(max_iter=10, random_state=0, sample_posterior=True)

results = []
for i in range(n_imputations):
    imputer.set_params(random_state=i)
    imputed = pd.DataFrame(imputer.fit_transform(df[numeric_cols]), columns=numeric_cols)
    for col in ['ARM', 'sex']:
        imputed[col] = df[col].values
    model = smf.logit(
        'outcome ~ C(ARM, Treatment(reference="Placebo")) + age', data=imputed
    ).fit(disp=0)
    results.append({'coef': model.params.iloc[1], 'se': model.bse.iloc[1]})

# Rubin's rules for pooling
pooled_coef = np.mean([r['coef'] for r in results])
within_var = np.mean([r['se']**2 for r in results])
between_var = np.var([r['coef'] for r in results], ddof=1)
total_var = within_var + (1 + 1 / n_imputations) * between_var
pooled_se = np.sqrt(total_var)
pooled_or = np.exp(pooled_coef)
pooled_ci = (np.exp(pooled_coef - 1.96 * pooled_se), np.exp(pooled_coef + 1.96 * pooled_se))
```

`sample_posterior=True` is essential. Without it, all m imputations produce nearly identical values, defeating the purpose of multiple imputation. This parameter draws from the posterior predictive distribution rather than using point estimates. Critical limitation: `sample_posterior=True` only works with `BayesianRidge` (the default estimator). If the estimator is changed (e.g., to `RandomForestRegressor`), the parameter is silently ignored and MI degenerates to single imputation.

Important: only impute covariates, not treatment assignment (which is fully determined by randomization) or the outcome (which creates circular dependency). Include outcome and treatment in the imputation model as predictors but exclude them from the set of imputed variables. For binary covariates, IterativeImputer treats them as continuous; consider `miceforest` for proper handling of mixed types.

The imputation model must be at least as flexible as the analysis model (congeniality). If the analysis includes treatment-by-covariate interactions, the imputation model should include them too. Uncongenial imputation biases estimates and invalidates variance pooling.

Note: `IterativeImputer` remains experimental in scikit-learn (requires the `enable_iterative_imputer` import) and its API may change without a standard deprecation cycle. For reproducibility-critical analyses, consider `miceforest` as a stable alternative.

The CI uses z=1.96 as a simplification. Proper Rubin's rules use a t-distribution with degrees of freedom: `df = (m-1) * (1 + W / ((1+1/m)*B))^2`. With m=20 the z-approximation is adequate; with fewer imputations, the t-distribution provides better coverage. The adequate number of imputations depends on the fraction of missing information (FMI): m >= 100 * FMI as a rule of thumb. With 40% missingness and FMI ~0.3, m=30 is needed.

When missing data exceeds 40% on key variables, results should be treated as hypothesis-generating unless MCAR can be demonstrated. Regulatory-standard sensitivity analyses include LOCF (Last Observation Carried Forward), BOCF (Baseline Observation Carried Forward), and tipping point analysis. While LOCF/BOCF are biased under MAR, regulators (especially FDA) still request them. Tipping point analysis asks how extreme the missing outcomes would need to be to overturn the primary conclusion and is the standard MNAR sensitivity approach.

### Rubin's Rules Summary

| Component | Formula |
|-----------|---------|
| Pooled estimate | Mean of m estimates |
| Within-imputation variance | Mean of m variance estimates |
| Between-imputation variance | Variance of m estimates |
| Total variance | W + (1 + 1/m) * B |
| Fraction of missing info | (1 + 1/m) * B / T |

## Multiplicity for Co-Primary Endpoints

| Method | Approach | Conservatism |
|--------|---------|-------------|
| Bonferroni | Divide alpha by number of endpoints | Most conservative |
| Hierarchical (gatekeeping) | Test in pre-specified order; proceed only if previous significant | Moderate |
| Hochberg step-up | Ordered p-values compared to alpha/(m-k+1) | Less conservative than Bonferroni |

Hierarchical testing requires a pre-specified ordering. The ordering should be based on clinical importance, not expected effect size.

## ICH E9(R1) Estimands Framework

An estimand defines precisely what is being estimated. Four attributes: population, variable/endpoint, intercurrent events strategy, population-level summary.

| Strategy | Approach | Example |
|----------|---------|---------|
| Treatment policy | Include all data regardless | ITT analysis |
| Composite | Incorporate event into endpoint | Death = non-responder |
| Hypothetical | Estimate as if event did not occur | Effect if no discontinuation |
| Principal stratum | Subpopulation who would not experience event | Completers regardless of arm |
| While on treatment | Data only while on assigned treatment | Per-protocol-like |

The estimand must be specified BEFORE choosing the statistical method. A common error is choosing a method (LOCF, MMRM) then retrofitting the estimand to match. The estimand drives the analysis, not vice versa.

## CONSORT Reporting Checklist

Key statistical requirements for trial reports:

- Effect estimates with confidence intervals and p-values for all primary and secondary endpoints
- Sample size calculation with stated assumptions (alpha, power, expected effect, dropout rate)
- Flow diagram showing participant numbers at each stage (screened, randomized, allocated, analyzed)
- Pre-specified methods for subgroup and adjusted analyses described before results
- Harms assessment methods specified, including how adverse events were collected and coded

```python
# Flow diagram counts
flow = {
    'screened': len(screening_log),
    'eligible': len(screening_log[screening_log['eligible']]),
    'randomized': len(dm),
    'allocated_drug': len(dm[dm['ARM'] == 'Drug']),
    'allocated_placebo': len(dm[dm['ARM'] == 'Placebo']),
    'completed_drug': len(dm[(dm['ARM'] == 'Drug') & dm['USUBJID'].isin(completers)]),
    'completed_placebo': len(dm[(dm['ARM'] == 'Placebo') & dm['USUBJID'].isin(completers)]),
    'analyzed_itt': len(itt),
}
for stage, count in flow.items():
    print(f'{stage}: {count}')
```

## Common Pitfalls

- **sample_posterior=False in IterativeImputer:** Underestimates uncertainty across imputations, producing artificially narrow confidence intervals.
- **Per-protocol as primary:** Inflates effect estimates and violates ICH E9 recommendations. ITT should be the primary analysis.
- **Missing confidence intervals:** Effect sizes without uncertainty measures are uninterpretable. Always report CI alongside point estimates.
- **Method before estimand:** The estimand framework (ICH E9 R1) requires defining what is estimated before selecting how to estimate it.
- **p-value tunnel vision:** Statistical significance alone does not establish clinical relevance. For each endpoint, state the pre-specified minimum clinically important difference (MCID) and whether the observed effect and its CI exceed it. An OR of 1.08 (95% CI 1.01-1.15, p=0.02) may be statistically significant but clinically negligible if the MCID was OR >= 1.25.
- **Post-hoc missing data handling:** The missing data strategy should be pre-specified in the SAP. Post-hoc approaches invite bias allegations.

## Related Skills

- clinical-biostatistics/cdisc-data-handling - Load and prepare the analysis dataset
- clinical-biostatistics/logistic-regression - Primary analysis models
- clinical-biostatistics/effect-measures - Compute and interpret treatment effects
- clinical-biostatistics/subgroup-analysis - Pre-specified subgroup analyses
- reporting/rmarkdown-reports - Generate formatted statistical reports
- experimental-design/multiple-testing - General multiplicity correction

---
name: bio-workflows-clinical-trial-pipeline
description: End-to-end clinical trial analysis workflow from CDISC data loading through statistical testing to regulatory-compliant reporting. Covers data preparation, logistic regression, categorical tests, subgroup analysis, and Table 1 generation. Use when performing a complete analysis of clinical trial data.
tool_type: python
primary_tool: statsmodels
workflow: true
depends_on:
  - clinical-biostatistics/cdisc-data-handling
  - clinical-biostatistics/logistic-regression
  - clinical-biostatistics/categorical-tests
  - clinical-biostatistics/effect-measures
  - clinical-biostatistics/subgroup-analysis
  - clinical-biostatistics/trial-reporting
qc_checkpoints:
  - after_data_prep: "One row per USUBJID, no duplicate subjects, treatment arms balanced"
  - after_primary_analysis: "Model converged, no separation warnings, pseudo-R2 reported"
  - after_subgroup: "Interaction tests run, multiplicity adjusted, forest plot generated"
  - after_reporting: "Table 1 complete, missing data documented, CONSORT items addressed"
---

## Version Compatibility

Reference examples tested with: statsmodels 0.14+, scipy 1.12+, tableone 0.9+, pyreadstat 1.2+, pandas 2.1+, numpy 1.26+, matplotlib 3.8+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Clinical Trial Analysis Pipeline

**"Analyze my clinical trial data end to end"** -> Load CDISC domain tables, prepare a subject-level analysis dataset, run primary statistical models, perform subgroup analyses, and generate regulatory-compliant tables and figures.

Complete workflow for clinical trial statistical analysis from raw data to publication-ready results.

### Scientific Reasoning Framework

Before executing any analysis step, establish the causal framework. For an RCT, randomization justifies causal interpretation of the primary analysis, but subgroup analyses and observational comparisons within the trial (e.g., adherence effects) do not inherit this protection. Key decisions requiring scientific judgment at each step: (1) data preparation -- which aggregation strategy matches the estimand, (2) covariate selection -- include confounders and prognostic factors from the SAP, exclude mediators and colliders, (3) subgroup analysis -- test only biologically motivated interactions, (4) missing data -- link DS domain reasons to the assumed mechanism before choosing a method. The workflow below provides the technical steps; the scientific reasoning at each decision point determines whether the results are valid.

## Workflow Overview

```
CDISC Domain Files (DM, AE, EX, LB)
    |
    v
[1. Data Preparation] ----> Subject-level dataset with outcomes and covariates
    |
    v
[2. Table 1] ------------> Baseline characteristics by treatment arm
    |
    v
[3. Primary Analysis] ---> Logistic regression with OR extraction
    |
    v
[4. Categorical Tests] --> Chi-square / Fisher's exact for key associations
    |
    v
[5. Subgroup Analysis] --> Interaction terms, stratified ORs, forest plot
    |
    v
[6. Missing Data] -------> Multiple imputation sensitivity analysis
    |
    v
Results tables and figures
```

## Step 1: Data Preparation

**Goal:** Create a single subject-level analysis dataset from CDISC domain tables.

**Approach:** Load domain files, aggregate event-level data to one row per subject, merge on USUBJID, and code the outcome variable.

```python
import pandas as pd
import pyreadstat

dm, _ = pyreadstat.read_xport('dm.xpt')
ae, _ = pyreadstat.read_xport('ae.xpt')

# Aggregate: did each subject have the target adverse event?
target_ae = ae[ae['AEDECOD'] == 'COVID-19']
severity_map = {'MILD': 1, 'MODERATE': 2, 'SEVERE': 3, 'LIFE THREATENING': 4, 'FATAL': 5}
target_ae['AESEV_NUM'] = target_ae['AESEV'].map(severity_map)
had_event = target_ae.groupby('USUBJID')['AESEV_NUM'].max().reset_index()
had_event.columns = ['USUBJID', 'EVENT_SEVERITY']

analysis = dm[['USUBJID', 'ARM', 'ARMCD', 'AGE', 'SEX']].merge(had_event, on='USUBJID', how='left')
analysis['HAD_EVENT'] = analysis['EVENT_SEVERITY'].notna().astype(int)
analysis['TREATMENT'] = (analysis['ARMCD'] != 'PLACEBO').astype(int)
```

**QC Checkpoint:** Verify one row per USUBJID, no unexpected duplicates, treatment arms are present and reasonably balanced.

```python
assert analysis['USUBJID'].is_unique, 'Duplicate subjects detected'
print(analysis['ARM'].value_counts())
```

## Step 2: Table 1 Baseline Characteristics

**Goal:** Summarize demographics and baseline variables by treatment arm.

**Approach:** Use TableOne to generate a publication-ready table with p-values and standardized mean differences.

```python
from tableone import TableOne

columns = ['AGE', 'SEX', 'RACE']
categorical = ['SEX', 'RACE']
table1 = TableOne(analysis, columns=columns, categorical=categorical,
                  groupby='ARM', pval=True, smd=True, missing=True)
print(table1.tabulate(tablefmt='github'))
```

Interpret SMD > 0.1 as meaningful imbalance rather than relying on p-values, which test whether randomization worked (a known mechanism, not a hypothesis).

## Step 3: Primary Analysis -- Logistic Regression

**Goal:** Estimate the treatment effect on the binary outcome as an adjusted odds ratio.

**Approach:** Fit a logistic regression with explicit reference category and clinically relevant covariates, then exponentiate coefficients to obtain ORs.

```python
import statsmodels.formula.api as smf
import numpy as np

model = smf.logit(
    'HAD_EVENT ~ C(ARM, Treatment(reference="Placebo")) + AGE + C(SEX)',
    data=analysis
).fit()

or_table = pd.DataFrame({
    'OR': np.exp(model.params),
    'Lower_CI': np.exp(model.conf_int()[0]),
    'Upper_CI': np.exp(model.conf_int()[1]),
    'p_value': model.pvalues
})
print(or_table)
print(f'McFadden pseudo-R2: {model.prsquared:.4f}')
```

**QC Checkpoint:** Verify model converged (no warnings), check for separation (coefficients > 10 or SE > 100), report pseudo-R-squared (McFadden > 0.2 is excellent; do not compare across pseudo-R2 types).

## Step 4: Categorical Tests

**Goal:** Test the crude association between treatment and outcome using contingency tables.

**Approach:** Build a 2x2 table, check expected cell counts, and choose chi-square or Fisher's exact accordingly.

```python
from scipy.stats import chi2_contingency, fisher_exact

table = pd.crosstab(analysis['ARM'], analysis['HAD_EVENT'])
chi2, p, dof, expected = chi2_contingency(table, correction=False)

if (expected < 5).any():
    _, p = fisher_exact(table.values)
    print(f'Fisher exact p = {p:.4f}')
else:
    print(f'Chi-square p = {p:.4f} (chi2 = {chi2:.2f}, dof = {dof})')
```

## Step 5: Subgroup Analysis

**Goal:** Test whether the treatment effect varies across pre-specified subgroups.

**Approach:** Fit a model with an interaction term, extract subgroup-specific ORs, adjust for multiplicity, and visualize with a forest plot.

```python
import matplotlib.pyplot as plt

# Interaction test
interaction_model = smf.logit(
    'HAD_EVENT ~ C(ARM, Treatment(reference="Placebo")) * C(SUBGROUP)',
    data=analysis
).fit()

# Subgroup-specific ORs
labels, ors, lowers, uppers = [], [], [], []
for group in analysis['SUBGROUP'].unique():
    sub = analysis[analysis['SUBGROUP'] == group]
    sub_model = smf.logit(
        'HAD_EVENT ~ C(ARM, Treatment(reference="Placebo"))',
        data=sub
    ).fit(disp=0)
    or_val = np.exp(sub_model.params.iloc[1])
    ci = np.exp(sub_model.conf_int().iloc[1])
    labels.append(group)
    ors.append(or_val)
    lowers.append(ci[0])
    uppers.append(ci[1])

# Multiplicity correction for subgroup p-values
from statsmodels.stats.multitest import multipletests
sub_pvals = [smf.logit('HAD_EVENT ~ C(ARM, Treatment(reference="Placebo"))',
             data=analysis[analysis['SUBGROUP'] == g]).fit(disp=0).pvalues.iloc[1]
             for g in labels]
_, adjusted_pvals, _, _ = multipletests(sub_pvals, method='holm')

# Forest plot
fig, ax = plt.subplots(figsize=(8, 5))
y_pos = range(len(labels))
ax.errorbar(ors, y_pos,
            xerr=[np.array(ors) - np.array(lowers), np.array(uppers) - np.array(ors)],
            fmt='D', color='black', capsize=3, markersize=5)
ax.axvline(x=1.0, color='gray', linestyle='--', linewidth=0.8)
ax.set_yticks(y_pos)
ax.set_yticklabels(labels)
ax.set_xlabel('Odds Ratio (95% CI)')
ax.set_xscale('log')
plt.tight_layout()
plt.savefig('forest_plot.png', dpi=150)
```

**QC Checkpoint:** Interaction p-value reported. Multiplicity correction applied if testing multiple subgroups. Forest plot shows overall estimate for context.

## Step 6: Missing Data Sensitivity Analysis

**Goal:** Assess robustness of the primary result under multiple imputation.

**Approach:** Impute missing covariates only (never treatment or outcome in an RCT), fit the primary model on each imputed dataset, and pool via Rubin's rules.

```python
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer

n_imputations = 20
covariate_cols = ['AGE']  # Only impute covariates, not treatment or outcome
mi_data = analysis.dropna(subset=['HAD_EVENT', 'TREATMENT']).copy()

results = []
for i in range(n_imputations):
    imputer = IterativeImputer(max_iter=10, random_state=i, sample_posterior=True)
    imputed_cov = pd.DataFrame(imputer.fit_transform(mi_data[covariate_cols]),
                               columns=covariate_cols, index=mi_data.index)
    imputed_cov['HAD_EVENT'] = mi_data['HAD_EVENT'].values
    imputed_cov['TREATMENT'] = mi_data['TREATMENT'].values
    model_imp = smf.logit('HAD_EVENT ~ TREATMENT + AGE', data=imputed_cov).fit(disp=0)
    results.append({'coef': model_imp.params['TREATMENT'], 'se': model_imp.bse['TREATMENT']})

pooled_coef = np.mean([r['coef'] for r in results])
within_var = np.mean([r['se']**2 for r in results])
between_var = np.var([r['coef'] for r in results], ddof=1)
total_var = within_var + (1 + 1/n_imputations) * between_var
pooled_or = np.exp(pooled_coef)
pooled_ci = (np.exp(pooled_coef - 1.96 * np.sqrt(total_var)),
             np.exp(pooled_coef + 1.96 * np.sqrt(total_var)))
print(f'Pooled OR: {pooled_or:.3f} ({pooled_ci[0]:.3f}-{pooled_ci[1]:.3f})')
```

**QC Checkpoint:** Compare pooled OR and CI with the complete-case primary analysis. Large discrepancies suggest missing data may not be MCAR. Document the comparison.

## Result Reporting Checklist

- [ ] Table 1 with baseline characteristics by arm (SMD > 0.1 flagged)
- [ ] Primary analysis OR with 95% CI and p-value
- [ ] Analysis population defined (ITT/PP/Safety)
- [ ] Missing data rate reported per variable
- [ ] Sensitivity analysis under multiple imputation
- [ ] Subgroup forest plot with interaction p-values
- [ ] Multiplicity adjustment method stated
- [ ] CONSORT flow diagram numbers available

## Related Skills

- clinical-biostatistics/cdisc-data-handling - Detailed CDISC domain handling and SUPPQUAL pivoting
- clinical-biostatistics/logistic-regression - Advanced regression options (ordinal, Firth's)
- clinical-biostatistics/categorical-tests - CMH stratified tests and McNemar's
- clinical-biostatistics/effect-measures - NNT/NNH and non-collapsibility
- clinical-biostatistics/subgroup-analysis - Breslow-Day homogeneity and RERI
- clinical-biostatistics/trial-reporting - Estimands framework and ICH E9 compliance

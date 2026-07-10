---
name: bio-workflows-clinical-trial-pipeline
description: End-to-end clinical trial analysis workflow from CDISC SDTM/ADaM loading through ICH E9(R1) estimand-driven primary analysis to CONSORT 2025 regulatory-compliant reporting. Covers data preparation, FDA 2023 marginal vs conditional logistic regression, categorical tests with Boschloo, modern HTE/subgroup methods, missing-data sensitivity (MMRM, reference-based MI, Permutt tipping point), graphical multiplicity (Bretz-Maurer), survival analysis (Cox/RMST/competing risks) when applicable, and Table 1. Use when performing a complete analysis of clinical trial data.
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
  - clinical-biostatistics/missing-data-sensitivity
  - clinical-biostatistics/multiplicity-graphical
  - clinical-biostatistics/survival-analysis
  - clinical-biostatistics/power-and-sample-size
qc_checkpoints:
  - after_estimand_definition: "ICH E9(R1) 5 attributes pre-specified in SAP: treatment, population, endpoint, summary measure, ICE handling strategy"
  - after_data_prep: "One row per USUBJID, no duplicate subjects, treatment arms balanced, DS domain tabulated for dropout patterns by arm"
  - after_primary_analysis: "Model converged, no separation warnings, marginal RD via g-computation reported as primary per FDA 2023; conditional OR as supportive"
  - after_subgroup: "Interaction tests run via single model (not per-subgroup p-comparisons), graphical multiplicity adjustment via gMCP, forest plot generated"
  - after_missing_data: "Per ICH E9(R1) ICE strategy: MMRM/MAR or reference-based MI (J2R/CR/CIR); Permutt tipping-point delta reported in residual SD units"
  - after_reporting: "Table 1 with SMD, missing data per CONSORT 2025 item 21c, harms per item 15, estimand statement per ICH E9(R1)"
---

## Version Compatibility

Reference examples tested with: statsmodels 0.14+, scipy 1.12+, tableone 0.9+, pyreadstat 1.2+, pandas 2.1+, numpy 1.26+, matplotlib 3.8+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Clinical Trial Analysis Pipeline

**"Analyze my clinical trial data end to end"** -> Load CDISC domain tables, prepare a subject-level analysis dataset, run primary statistical models, perform subgroup analyses, and generate regulatory-compliant tables and figures.

This is a workflow skill: it owns the chaining decisions and hand-offs, not the internals of any one step.

## The governing principle

A clinical-trial analysis is a commit-then-execute pipeline: its trustworthiness is decided by whether the analysis was locked BEFORE the data was seen, and every seam failure is a wrong-but-silent handoff that answers a different question with no error thrown.

1. **The ESTIMAND, locked in the SAP before unblinding, is the made-once commitment everything inherits.** ICH E9(R1) defines it by 5 attributes: treatment condition, population, endpoint, intercurrent-event (ICE) handling strategy, and population-level summary. It is the precise definition of "what treatment effect the trial estimates", committed before analysis so the estimator is matched to it rather than chosen to flatter the data.
2. **The estimator is a CONSEQUENCE of the estimand's ICE strategy, not a free modeling choice.** A treatment-policy strategy needs all post-ICE data (MMRM/treatment-policy, or reference-based MI if truly missing); a hypothetical strategy censors/models the post-ICE data as if the ICE had not occurred (MMRM under MAR); a composite strategy folds the ICE into the endpoint. Choosing MMRM vs reference-based MI vs g-computation to flatter the data is an estimand-estimator mismatch.
3. **The analysis population and all subgroups are defined a priori; the data never picks them.** Randomization licenses causal interpretation for the ITT/FAS primary only — PP and subgroups do NOT inherit that protection. Pre-specify every subgroup and the multiplicity graph in the SAP; post-hoc data-driven subgroups are hypothesis-generating only. This is the clinical analog of "don't call hits before CN correction".
4. **Baseline balance is not tested with p-values, and marginal is not conditional.** In a randomized trial any imbalance is by definition chance, and baseline hypothesis testing is incoherent (Senn 1994) — report SMD (>0.1 notable, a convention Senn does not state) and adjust via pre-specified ANCOVA, don't test-then-decide. For a binary endpoint the primary is the MARGINAL risk difference via g-computation (FDA 2023); a conditional OR is a different parameter (non-collapsibility), reported as supportive.

## Made-once commitments

| Commitment | Consequence inherited downstream |
|------------|----------------------------------|
| The estimand (5 ICH E9(R1) attributes, SAP-locked before unblinding) | Which question the trial answers; the estimator must match its ICE strategy |
| The estimator matched to the ICE strategy | Whether the analysis answers the locked question; a mismatch is silent |
| Analysis population set (ITT/FAS primary, PP sensitivity, Safety as-treated) a priori | Causal validity; ITT is randomization-protected, PP and subgroups are not |
| SDTM -> ADaM derivation (pre-specified, one-directional; ADTTE CNSR convention) | Every downstream model; CDISC CNSR=0 means EVENT, opposite of R/Python survival packages |

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
[3. Primary Analysis] ---> Marginal RD via g-computation (conditional OR supportive)
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
target_ae = ae[ae['AEDECOD'] == 'COVID-19'].copy()
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

**Approach:** Use TableOne to generate a baseline table by arm with standardized mean differences and explicit missingness. Omit the baseline p-value column: in a randomized trial any imbalance is by definition due to chance, so a baseline p-value tests a null already known to be true (Senn 1994; CONSORT 2010/2025). Report SMD for balance instead. Table-construction and export mechanics (gtsummary/tableone, Word export, gene-symbol-safe supplements) live in reporting/publication-tables.

```python
from tableone import TableOne

columns = ['AGE', 'SEX', 'RACE']
categorical = ['SEX', 'RACE']
table1 = TableOne(analysis, columns=columns, categorical=categorical,
                  groupby='ARM', pval=False, smd=True, missing=True)
print(table1.tabulate(tablefmt='github'))
```

Interpret SMD > 0.1 as meaningful imbalance. The response to a worrying imbalance on a prognostic covariate is to adjust for it (a pre-specified ANCOVA/model covariate), not to test it.

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

The logistic fit above yields the CONDITIONAL (adjusted) OR. When the estimand's summary measure is a marginal risk difference (the FDA 2023 primary for a binary endpoint), do NOT report this conditional OR as the primary effect. Compute the MARGINAL risk difference by g-computation: fit the covariate-adjusted model, predict each subject's outcome probability under both arms, average within arm, and difference; bootstrap the whole fit-and-predict for the CI. `examples/clinical_trial_pipeline.py` implements this; see clinical-biostatistics/logistic-regression for the estimator's assumptions and variance options. Report the conditional OR as supportive. The two differ by non-collapsibility and answer different questions.

**QC Checkpoint:** Verify model converged (no warnings), check for separation (coefficients > 10 or SE > 100), report pseudo-R-squared (McFadden > 0.2 is excellent; do not compare across pseudo-R2 types). Confirm the reported PRIMARY effect matches the estimand's summary measure (marginal RD via g-computation for a marginal estimand), not whichever the model emits by default.

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

**Approach:** Fit a model with an interaction term and test it with a single interaction LR test. Subgroup-specific ORs are DESCRIPTIVE ONLY -- do not correct their individual p-values, do not interpret them; the alpha allocated to the subgroup family is handled by the pre-specified gMCP graph, not by a post-hoc correction. Visualize with a forest plot.

```python
import matplotlib.pyplot as plt

# The HTE test is the treatment-by-subgroup INTERACTION, not per-subgroup significance.
# Fit the main-effects (restricted) and interaction (full) models, then LR-test the interaction.
main_model = smf.logit(
    'HAD_EVENT ~ C(ARM, Treatment(reference="Placebo")) + C(SUBGROUP)',
    data=analysis
).fit(disp=0)
interaction_model = smf.logit(
    'HAD_EVENT ~ C(ARM, Treatment(reference="Placebo")) * C(SUBGROUP)',
    data=analysis
).fit(disp=0)
# compare_lr_test exists only on linear-model results, NOT LogitResults -- compute the LR test by hand.
from scipy.stats import chi2
lr_stat = 2 * (interaction_model.llf - main_model.llf)
df_diff = int(interaction_model.df_model - main_model.df_model)
interaction_pval = chi2.sf(lr_stat, df_diff)
print(f'Treatment-by-subgroup interaction: LR chi2={lr_stat:.2f}, df={df_diff}, p={interaction_pval:.3f}')

# Subgroup-specific ORs are DESCRIPTIVE ONLY (to draw the forest plot). Do NOT interpret their
# individual p-values as evidence of a subgroup effect -- that is the invalid per-subgroup
# significance pattern; the interaction p-value above is the only valid HTE test.
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

**QC Checkpoint:** Interaction p-value reported from a single LR test. Alpha for the subgroup family allocated in the pre-specified gMCP graph (not a post-hoc p-value correction on the descriptive subgroup ORs). Forest plot shows overall estimate for context.

## Step 6: Missing Data Sensitivity Analysis (per ICH E9(R1) and clinical-biostatistics/missing-data-sensitivity)

**Goal:** Assess robustness of the primary result under the pre-specified ICE strategy with both MAR primary and MNAR sensitivity analyses.

**Approach:** First examine DS (Disposition) domain for differential dropout patterns; if dropout differs by arm, MAR is suspect and reference-based MI is required as primary. Otherwise, fit MMRM under MAR with Rubin's-rules pooling for continuous endpoints, or g-computation with bootstrap for binary. Always run Permutt 2016 tipping-point sensitivity in residual SD units.

```python
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer

n_imputations = 20   # practical starting count; von Hippel 2020 shows required m scales with the fraction of missing information (two-stage rule), so raise it when FMI is high
# Impute AGE jointly WITH its predictors (SEX, TREATMENT, HAD_EVENT). A single-column imputer has
# no predictors, so sample_posterior draws are identical -> between-imputation variance = 0 and the
# MI collapses to complete-case. Including correlated columns makes the posterior draws actually vary.
impute_cols = ['AGE', 'TREATMENT', 'HAD_EVENT']   # numeric only; SEX ('M'/'F') would break IterativeImputer and is restored below
mi_data = analysis.dropna(subset=['HAD_EVENT', 'TREATMENT']).copy()

results = []
for i in range(n_imputations):
    imputer = IterativeImputer(max_iter=10, random_state=i, sample_posterior=True)
    imputed_cov = pd.DataFrame(imputer.fit_transform(mi_data[impute_cols]),
                               columns=impute_cols, index=mi_data.index)
    # Only AGE had missings; restore the observed discrete columns so they stay integer-valued.
    imputed_cov['HAD_EVENT'] = mi_data['HAD_EVENT'].values
    imputed_cov['TREATMENT'] = mi_data['TREATMENT'].values
    imputed_cov['SEX'] = mi_data['SEX'].values
    # Mirror the primary ADJUSTED model's RHS (AGE + SEX) so the pooled OR is comparable to it.
    model_imp = smf.logit('HAD_EVENT ~ TREATMENT + AGE + C(SEX)', data=imputed_cov).fit(disp=0)
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

## Result Reporting Checklist (CONSORT 2025 + ICH E9(R1) aligned)

- [ ] ICH E9(R1) estimand statement with 5 attributes pre-specified in SAP
- [ ] Table 1 with baseline characteristics by arm (SMD > 0.1 flagged; NOT p-values)
- [ ] Primary analysis: marginal RD via g-computation per FDA 2023 (binary) OR MMRM-MAR with Kenward-Roger (continuous longitudinal)
- [ ] Conditional OR/HR as supportive (different parameter than marginal due to non-collapsibility)
- [ ] Analysis populations defined: ITT (primary), FAS (with explicit exclusion criteria), PP (sensitivity), Safety (AE)
- [ ] Missing data per CONSORT 2025 item 21c: mechanism assumption, primary method, MNAR sensitivity (J2R/CR/CIR per Carpenter-Roger 2013)
- [ ] Permutt tipping-point delta reported in residual SD units
- [ ] Subgroup forest plot with INTERACTION p-values (not per-subgroup p-comparison); graphical multiplicity via gMCP
- [ ] Multiplicity adjustment method stated (CONSORT 2025 item 20; FDA Multiple Endpoints Final Oct 2022)
- [ ] CONSORT flow diagram numbers available
- [ ] Harms per CONSORT 2025 item 15 (absorbs CONSORT-Harms 2022)

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| The analysis answers a different question | Estimand-estimator mismatch (e.g. a hypothetical estimand analyzed treatment-policy) | Derive the estimator FROM the ICE strategy; MMRM/MI/g-computation are consequences of attribute (iv), not free choices |
| Spurious subgroup effect | Subgroups/sensitivity chosen after seeing the data | Pre-specify all subgroups + the multiplicity graph in the SAP; post-hoc is hypothesis-generating only |
| Effect loses causal validity | PP swapped in as primary because it "looks cleaner" | ITT/FAS primary (randomization-protected); PP is sensitivity; define both a priori |
| Survival results inverted | ADTTE CNSR sign flip (CDISC CNSR=0 means EVENT) | Convert the censoring indicator before passing to lifelines/survival |
| False "imbalance" conclusions | Baseline characteristics tested with p-values | Report SMD (>0.1 notable); adjust prognostic imbalance via pre-specified ANCOVA |
| Primary effect is the wrong parameter | Conditional OR reported as the marginal effect (non-collapsibility) | Marginal risk difference via g-computation is primary (FDA 2023); conditional OR supportive |
| Over-conservative inference (CIs too wide, power lost) | Stratification/randomization factors omitted from the analysis model | Include the stratification factors; omitting them biases SEs upward and Type-I below nominal (Kahan-Morris 2012) |

## When to Add Specialized Skills

This pipeline covers the typical binary-endpoint RCT workflow. For specific designs, add the corresponding specialized skill:

- **Time-to-event primary endpoint** (OS, PFS, DOR): add clinical-biostatistics/survival-analysis for Cox PH diagnostics, RMST under non-PH, competing risks via Fine-Gray vs cause-specific Cox, MaxCombo for delayed effects, informative censoring handling
- **Continuous longitudinal endpoint** (HbA1c at 24 weeks): clinical-biostatistics/missing-data-sensitivity for MMRM with Kenward-Roger via R mmrm; reference-based MI via R rbmi for MNAR sensitivity
- **Multiple primary or key secondary endpoints**: clinical-biostatistics/multiplicity-graphical for Bretz-Maurer graphical procedures via gMCP
- **Trial design / sample-size justification**: clinical-biostatistics/power-and-sample-size for Schoenfeld events, Lakatos under non-PH, FDA 2016 NI double discount, TOST equivalence
- **Adaptive trial** (group-sequential, SSR, platform): clinical-biostatistics/adaptive-designs for rpact/gsDesign, Mehta-Pocock promising zone, ICH E20 considerations
- **Bayesian primary inference or RWE comparator**: clinical-biostatistics/bayesian-trials for BOIN dose-finding, robust MAP priors via RBesT, EXNEX basket trials, psborrow2 for external controls

## Related Skills

- clinical-biostatistics/cdisc-data-handling - CDISC SDTM/ADaM, Pinnacle 21, Dataset-JSON, ADTTE CNSR conventions
- clinical-biostatistics/logistic-regression - FDA 2023 marginal vs conditional, g-computation, Brant test, Firth, Hauck-Donner
- clinical-biostatistics/categorical-tests - Boschloo, mid-p McNemar, Wilson/Newcombe/Miettinen-Nurminen CIs
- clinical-biostatistics/effect-measures - NNT Bender 2002 convention, profile likelihood, modified Poisson for RR
- clinical-biostatistics/subgroup-analysis - Causal forests, STEPP, SIDES, EXNEX, Yadlowsky RATE, EMA 2019 subgroup guideline
- clinical-biostatistics/trial-reporting - ICH E9(R1) 5 estimand strategies, Cro/Bartlett variance debate, CONSORT 2025
- clinical-biostatistics/missing-data-sensitivity - MMRM/Kenward-Roger, reference-based MI, Permutt tipping point
- clinical-biostatistics/multiplicity-graphical - Bretz-Maurer graphs, Goeman closed-testing admissibility
- clinical-biostatistics/survival-analysis - Cox/RMST/Fine-Gray/MaxCombo/recurrent events/interval censoring
- clinical-biostatistics/power-and-sample-size - Schoenfeld/Lakatos, NI double discount, crossover, MCID
- clinical-biostatistics/adaptive-designs - Group-sequential, SSR, RAR consensus, BOIN, platform trials
- clinical-biostatistics/bayesian-trials - MAP/EXNEX/RWE, FDA Bayesian Jan 2026 draft, psborrow2
- reporting/publication-tables - Table 1 construction (SMD not baseline p-values, show missingness) and Word/LaTeX export

## References

- ICH E9(R1) (2019) Addendum on Estimands and Sensitivity Analysis in Clinical Trials to the Guideline on Statistical Principles for Clinical Trials. International Council for Harmonisation. (the estimand's 5 attributes and 5 ICE strategies.)
- Kahan BC, Cro S, Li F, et al (2023) Eliminating ambiguous treatment effects using estimands. *American Journal of Epidemiology* 192:987-994. DOI 10.1093/aje/kwad036. (98% of trials do not articulate the estimand.)
- Senn S (1994) Testing for baseline balance in clinical trials. *Statistics in Medicine* 13:1715-1726. DOI 10.1002/sim.4780131703. (baseline SMD, not p-values.)
- Kahan BC, Morris TP (2012) Improper analysis of trials randomised using stratified blocks or minimisation. *Statistics in Medicine* 31:328-340. DOI 10.1002/sim.4431. (omitting stratification factors makes inference conservative: SEs biased upward, Type-I below nominal, power lost.)
- Carpenter JR, Roger JH, Kenward MG (2013) Analysis of longitudinal trials with protocol deviation: a framework for relevant, accessible assumptions, and inference via multiple imputation. *Journal of Biopharmaceutical Statistics* 23:1352-1371. DOI 10.1080/10543406.2013.834911. (reference-based MI.)

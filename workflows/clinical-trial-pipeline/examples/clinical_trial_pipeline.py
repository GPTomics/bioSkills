'''End-to-end clinical trial analysis pipeline'''
# Reference: statsmodels 0.14+, scipy 1.12+, tableone 0.9+, pandas 2.1+, numpy 1.26+ | Verify API if version differs

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import chi2_contingency, fisher_exact
from tableone import TableOne
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer

# --- Step 1: Data Preparation ---
dm = pd.read_csv('DM.csv')
ae = pd.read_csv('AE.csv')

target_ae = ae[ae['AEDECOD'] == 'COVID-19'].copy()
severity_map = {'MILD': 1, 'MODERATE': 2, 'SEVERE': 3, 'LIFE THREATENING': 4, 'FATAL': 5}
target_ae['AESEV_NUM'] = target_ae['AESEV'].map(severity_map)
had_event = target_ae.groupby('USUBJID')['AESEV_NUM'].max().reset_index()
had_event.columns = ['USUBJID', 'EVENT_SEVERITY']

analysis = dm[['USUBJID', 'ARM', 'ARMCD', 'AGE', 'SEX']].merge(had_event, on='USUBJID', how='left')
analysis['HAD_EVENT'] = analysis['EVENT_SEVERITY'].notna().astype(int)
analysis['TREATMENT'] = (analysis['ARMCD'] != 'PLACEBO').astype(int)

assert analysis['USUBJID'].is_unique, 'Duplicate subjects detected'
print(f'Subjects: {len(analysis)}, Events: {analysis["HAD_EVENT"].sum()}')
print(analysis['ARM'].value_counts())

# --- Step 2: Table 1 ---
columns = ['AGE', 'SEX']
table1 = TableOne(analysis, columns=columns, categorical=['SEX'],
                  groupby='ARM', pval=False, smd=True, missing=True)   # SMD, NOT baseline p-values (Senn 1994)
print(table1.tabulate(tablefmt='github'))

# --- Step 3: Primary Analysis ---
model = smf.logit('HAD_EVENT ~ C(ARM, Treatment(reference="Placebo")) + AGE + C(SEX)', data=analysis).fit()
or_table = pd.DataFrame({
    'OR': np.exp(model.params),
    'Lower_CI': np.exp(model.conf_int()[0]),
    'Upper_CI': np.exp(model.conf_int()[1]),
    'p_value': model.pvalues
})
print('\nSupportive Analysis - Conditional Odds Ratios:')
print(or_table.to_string(float_format='%.4f'))
print(f'McFadden pseudo-R2: {model.prsquared:.4f}')

# The conditional OR above is SUPPORTIVE. The FDA 2023 primary for a binary endpoint is the MARGINAL
# risk difference: the OR is non-collapsible, so it answers a different question than the RD. Estimate
# by g-computation -- predict every subject's risk under BOTH arms, average within arm, difference.
REFERENCE_ARM = 'Placebo'
ACTIVE_ARM = next(a for a in analysis['ARM'].unique() if a != REFERENCE_ARM)

def marginal_rd(fitted, data):
    return fitted.predict(data.assign(ARM=ACTIVE_ARM)).mean() - fitted.predict(data.assign(ARM=REFERENCE_ARM)).mean()

rd = marginal_rd(model, analysis)

# Percentile bootstrap: refit on each resample so the CI reflects estimation of the outcome model too.
# 1000 resamples is the conventional floor for a stable 2.5/97.5 percentile.
rng = np.random.default_rng(42)
boot = []
for _ in range(1000):
    # reset_index: sampling with replacement duplicates index labels, and whenever a covariate is
    # missing patsy drops those rows, forcing predict() to reindex -- which raises on duplicate labels.
    resample = analysis.iloc[rng.integers(0, len(analysis), len(analysis))].reset_index(drop=True)
    refit = smf.logit('HAD_EVENT ~ C(ARM, Treatment(reference="Placebo")) + AGE + C(SEX)', data=resample).fit(disp=0)
    boot.append(marginal_rd(refit, resample))
lo, hi = np.percentile(boot, [2.5, 97.5])
print(f'\nPrimary Analysis - marginal risk difference ({ACTIVE_ARM} vs {REFERENCE_ARM}), g-computation:')
print(f'RD = {rd:+.4f}  95% bootstrap CI [{lo:+.4f}, {hi:+.4f}]')

# --- Step 4: Chi-square Test ---
ct = pd.crosstab(analysis['ARM'], analysis['HAD_EVENT'])
chi2, p, dof, expected = chi2_contingency(ct, correction=False)
if (expected < 5).any():
    _, p = fisher_exact(ct.values)
    print(f'\nFisher exact p = {p:.4f}')
else:
    print(f'\nChi-square: chi2={chi2:.2f}, p={p:.4f}, dof={dof}')

# --- Step 5: Missing Data Sensitivity ---
n_imputations = 20   # practical starting count; von Hippel 2020 shows required m scales with the fraction of missing information (two-stage rule), so raise it when FMI is high
# Impute AGE jointly WITH its predictors (SEX, TREATMENT, HAD_EVENT). A single-column imputer has no
# predictors, so sample_posterior draws are identical -> between-imputation variance 0 and MI collapses
# to complete-case. Correlated columns make the posterior draws vary.
impute_cols = ['AGE', 'TREATMENT', 'HAD_EVENT']   # numeric only; SEX ('M'/'F') would break IterativeImputer and is restored below
mi_analysis = analysis.dropna(subset=['HAD_EVENT', 'TREATMENT']).copy()
mi_results = []
for i in range(n_imputations):
    imputer = IterativeImputer(max_iter=10, random_state=i, sample_posterior=True)
    imputed_covariates = pd.DataFrame(
        imputer.fit_transform(mi_analysis[impute_cols]),
        columns=impute_cols, index=mi_analysis.index
    )
    imputed = imputed_covariates.copy()
    imputed['HAD_EVENT'] = mi_analysis['HAD_EVENT'].values
    imputed['TREATMENT'] = mi_analysis['TREATMENT'].values
    imputed['SEX'] = mi_analysis['SEX'].values
    # Mirror the primary ADJUSTED model's RHS (AGE + SEX) so the pooled OR is comparable to it.
    mi_model = smf.logit('HAD_EVENT ~ TREATMENT + AGE + C(SEX)', data=imputed).fit(disp=0)
    mi_results.append({'coef': mi_model.params['TREATMENT'], 'se': mi_model.bse['TREATMENT']})

pooled_coef = np.mean([r['coef'] for r in mi_results])
within_var = np.mean([r['se']**2 for r in mi_results])
between_var = np.var([r['coef'] for r in mi_results], ddof=1)
total_var = within_var + (1 + 1/n_imputations) * between_var
pooled_or = np.exp(pooled_coef)
pooled_ci = (np.exp(pooled_coef - 1.96 * np.sqrt(total_var)), np.exp(pooled_coef + 1.96 * np.sqrt(total_var)))
print(f'\nMultiple Imputation: OR={pooled_or:.3f} ({pooled_ci[0]:.3f}-{pooled_ci[1]:.3f})')

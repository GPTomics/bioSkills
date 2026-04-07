'''Generate Table 1, define analysis populations, and demonstrate multiple imputation for a clinical trial'''
# Reference: tableone 0.9+, statsmodels 0.14+, sklearn 1.4+, pandas 2.1+ | Verify API if version differs

import numpy as np
import pandas as pd
from tableone import TableOne
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
import statsmodels.formula.api as smf

np.random.seed(42)
n = 400
df = pd.DataFrame({
    'USUBJID': [f'STUDY01-{i:04d}' for i in range(n)],
    'ARM': np.random.choice(['Drug', 'Placebo'], n),
    'age': np.random.normal(58, 12, n).round(1),
    'sex': np.random.choice(['M', 'F'], n, p=[0.55, 0.45]),
    'race': np.random.choice(['White', 'Black', 'Asian', 'Other'], n, p=[0.6, 0.2, 0.15, 0.05]),
    'bmi': np.random.normal(27, 5, n).round(1),
    'baseline_score': np.random.normal(50, 10, n).round(1),
    'disease_stage': np.random.choice(['Mild', 'Moderate', 'Severe'], n, p=[0.3, 0.5, 0.2]),
})

base_prob = 0.35
df['outcome'] = np.random.binomial(1, base_prob + (df['ARM'] == 'Drug').astype(float) * -0.12)
df['received_dose'] = np.random.binomial(1, 0.95, n)
df['completed'] = np.random.binomial(1, 0.85, n)
df['protocol_violation'] = np.random.binomial(1, 0.08, n)

missing_mask = np.random.rand(n) < 0.12
df.loc[missing_mask, 'baseline_score'] = np.nan
missing_outcome = np.random.rand(n) < 0.08
df.loc[missing_outcome, 'outcome'] = np.nan

print('=== Table 1: Baseline Characteristics ===')
columns = ['age', 'sex', 'race', 'bmi', 'baseline_score', 'disease_stage']
categorical = ['sex', 'race', 'disease_stage']
table1 = TableOne(df, columns=columns, categorical=categorical, groupby='ARM', pval=True, smd=True, missing=True, overall=True)
print(table1.tabulate(tablefmt='github'))

print('\n=== Analysis Populations ===')
itt = df.copy()
pp = df[(df['completed'] == 1) & (df['protocol_violation'] == 0)]
safety = df[df['received_dose'] == 1]

for name, pop in [('ITT', itt), ('Per-Protocol', pp), ('Safety', safety)]:
    arm_counts = pop['ARM'].value_counts().to_dict()
    print(f'{name}: n={len(pop)}, Drug={arm_counts.get("Drug", 0)}, Placebo={arm_counts.get("Placebo", 0)}')

print(f'\nMissing data: outcome={df["outcome"].isna().sum()} ({df["outcome"].isna().mean()*100:.1f}%), baseline_score={df["baseline_score"].isna().sum()} ({df["baseline_score"].isna().mean()*100:.1f}%)')

print('\n=== Multiple Imputation (Rubin\'s Rules) ===')
analysis_df = itt.dropna(subset=['ARM', 'outcome']).copy()
covariate_cols = ['age', 'bmi', 'baseline_score']
n_imputations = 20

results = []
for i in range(n_imputations):
    imputer = IterativeImputer(max_iter=10, random_state=i, sample_posterior=True)
    imputed_covariates = pd.DataFrame(
        imputer.fit_transform(analysis_df[covariate_cols]),
        columns=covariate_cols, index=analysis_df.index
    )
    imputed = imputed_covariates.copy()
    imputed['outcome'] = analysis_df['outcome'].values
    imputed['ARM'] = analysis_df['ARM'].values
    model = smf.logit('outcome ~ C(ARM, Treatment(reference="Placebo")) + age + bmi', data=imputed).fit(disp=0)
    results.append({'coef': model.params.iloc[1], 'se': model.bse.iloc[1]})

pooled_coef = np.mean([r['coef'] for r in results])
within_var = np.mean([r['se']**2 for r in results])
between_var = np.var([r['coef'] for r in results], ddof=1)
total_var = within_var + (1 + 1 / n_imputations) * between_var
pooled_se = np.sqrt(total_var)
pooled_or = np.exp(pooled_coef)
pooled_ci_lower = np.exp(pooled_coef - 1.96 * pooled_se)
pooled_ci_upper = np.exp(pooled_coef + 1.96 * pooled_se)
fmi = (1 + 1 / n_imputations) * between_var / total_var

print(f'Pooled OR (Drug vs Placebo): {pooled_or:.3f} ({pooled_ci_lower:.3f}-{pooled_ci_upper:.3f})')
print(f'Pooled coefficient: {pooled_coef:.4f} (SE={pooled_se:.4f})')
print(f'Fraction of missing information: {fmi:.3f}')

print('\n=== Complete-Case Sensitivity Analysis ===')
cc = analysis_df.dropna(subset=covariate_cols)
cc_model = smf.logit('outcome ~ C(ARM, Treatment(reference="Placebo")) + age + bmi', data=cc).fit(disp=0)
cc_or = np.exp(cc_model.params.iloc[1])
cc_ci = np.exp(cc_model.conf_int().iloc[1])
print(f'Complete-case OR: {cc_or:.3f} ({cc_ci[0]:.3f}-{cc_ci[1]:.3f}), n={len(cc)}')
print(f'MI OR:            {pooled_or:.3f} ({pooled_ci_lower:.3f}-{pooled_ci_upper:.3f}), n={len(analysis_df)}')

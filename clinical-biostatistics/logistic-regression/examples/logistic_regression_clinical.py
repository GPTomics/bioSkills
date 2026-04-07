'''Logistic regression for a binary clinical trial endpoint with covariate adjustment'''
# Reference: statsmodels 0.14+, pandas 2.1+ | Verify API if version differs

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

df = pd.read_csv('analysis_dataset.csv')
print(f'Loaded {len(df)} subjects')
print(f'Treatment arms: {df["ARM"].value_counts().to_dict()}')
print(f'Outcome prevalence: {df["outcome"].mean():.1%}')

print(f'\nOutcome by arm:')
cross = pd.crosstab(df['ARM'], df['outcome'], margins=True)
print(cross)

model = smf.logit(
    'outcome ~ C(ARM, Treatment(reference="Placebo")) + age + C(sex)',
    data=df
).fit()

print(f'\n=== Model Summary ===')
print(f'Pseudo R-squared (McFadden): {model.prsquared:.3f}')
print(f'Log-likelihood: {model.llf:.1f}')
print(f'AIC: {model.aic:.1f}')

or_table = pd.DataFrame({
    'OR': np.exp(model.params),
    'Lower_CI': np.exp(model.conf_int()[0]),
    'Upper_CI': np.exp(model.conf_int()[1]),
    'p_value': model.pvalues
})
or_table = or_table.drop('Intercept', errors='ignore')

print(f'\n=== Odds Ratios ===')
for idx, row in or_table.iterrows():
    sig = '*' if row['p_value'] < 0.05 else ''
    print(f'{idx}: OR={row["OR"]:.3f} (95% CI: {row["Lower_CI"]:.3f}-{row["Upper_CI"]:.3f}), p={row["p_value"]:.4f} {sig}')

or_table.to_csv('odds_ratios.csv')
print(f'\nSaved odds_ratios.csv')

print(f'\n=== Prediction Table (rows=actual, cols=predicted) ===')
print(model.pred_table())

y_pred = model.predict()

from sklearn.metrics import roc_auc_score
auc = roc_auc_score(df['outcome'], y_pred)
print(f'\nROC-AUC: {auc:.3f}')

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

hl_stat, hl_p = hosmer_lemeshow(df['outcome'].values, y_pred.values)
print(f'Hosmer-Lemeshow: chi2={hl_stat:.2f}, p={hl_p:.4f}')
print('(H-L is supplementary; use calibration plots as primary calibration assessment)')

max_coef = model.params.abs().max()
max_se = model.bse.max()
if max_coef > 10 or max_se > 100:
    print(f'\nWARNING: Large coefficient ({max_coef:.1f}) or SE ({max_se:.1f}) detected.')
    print('Possible separation -- consider Firth penalized regression.')

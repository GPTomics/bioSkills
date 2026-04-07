'''Compute treatment effect measures from a clinical trial dataset'''
# Reference: statsmodels 0.14+, numpy 1.26+, pandas 2.1+, matplotlib 3.8+ | Verify API if version differs

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.contingency_tables import Table2x2
import matplotlib.pyplot as plt

np.random.seed(42)
n = 300
df = pd.DataFrame({
    'ARM': np.random.choice(['Active', 'Placebo'], n),
    'age': np.random.normal(55, 12, n).astype(int),
    'sex': np.random.choice(['M', 'F'], n)
})
logit_p = -1.5 + 0.8 * (df['ARM'] == 'Active').astype(int) + 0.02 * df['age']
df['outcome'] = (np.random.uniform(size=n) < 1 / (1 + np.exp(-logit_p))).astype(int)

print(f'Loaded {len(df)} subjects')
print(f'Treatment arms: {df["ARM"].value_counts().to_dict()}')
print(f'Outcome prevalence: {df["outcome"].mean():.1%}')

cross = pd.crosstab(df['ARM'], df['outcome'])
cross = cross[[1, 0]]
print(f'\nContingency table (outcome=1 first for Table2x2):\n{cross}')

table_2x2 = cross.values
t = Table2x2(table_2x2)
print(f'\n=== Crude Measures (Table2x2) ===')
print(f'Odds Ratio: {t.oddsratio:.3f}')
or_ci = t.oddsratio_confint()
print(f'OR 95% CI: ({or_ci[0]:.3f}, {or_ci[1]:.3f})')
print(f'Log OR SE: {t.log_oddsratio_se:.4f}')
print(f'Risk Ratio: {t.riskratio:.3f}')
rr_ci = t.riskratio_confint()
print(f'RR 95% CI: ({rr_ci[0]:.3f}, {rr_ci[1]:.3f})')

model = smf.logit('outcome ~ C(ARM, Treatment(reference="Placebo")) + age + C(sex)', data=df).fit()
or_table = pd.DataFrame({
    'OR': np.exp(model.params),
    'Lower_CI': np.exp(model.conf_int()[0]),
    'Upper_CI': np.exp(model.conf_int()[1]),
    'p_value': model.pvalues
})
or_table = or_table.drop('Intercept', errors='ignore')
print(f'\n=== Adjusted Odds Ratios (Logistic Regression) ===')
for idx, row in or_table.iterrows():
    sig = '*' if row['p_value'] < 0.05 else ''
    print(f'{idx}: OR={row["OR"]:.3f} ({row["Lower_CI"]:.3f}-{row["Upper_CI"]:.3f}), p={row["p_value"]:.4f} {sig}')

def nnt_from_or(odds_ratio, baseline_risk):
    baseline_odds = baseline_risk / (1 - baseline_risk)
    treatment_odds = baseline_odds * odds_ratio
    treatment_risk = treatment_odds / (1 + treatment_odds)
    arr = abs(baseline_risk - treatment_risk)
    return int(np.ceil(1 / arr)) if arr > 0 else float('inf')

crude_or = t.oddsratio
baseline_risk = df[df['ARM'] == 'Placebo']['outcome'].mean()
nnt = nnt_from_or(crude_or, baseline_risk)
print(f'\n=== Number Needed to Treat ===')
print(f'Baseline risk (Placebo): {baseline_risk:.1%}')
print(f'Crude OR: {crude_or:.3f}')
print(f'NNT: {nnt}')

subgroups = {'Male': df[df['sex'] == 'M'], 'Female': df[df['sex'] == 'F'],
             'Age < 55': df[df['age'] < 55], 'Age >= 55': df[df['age'] >= 55],
             'Overall': df}
labels, ors, lower_cis, upper_cis = [], [], [], []
for label, sub_df in subgroups.items():
    sub_cross = pd.crosstab(sub_df['ARM'], sub_df['outcome'])
    sub_cross = sub_cross[[1, 0]].values if set([0, 1]).issubset(sub_cross.columns) else sub_cross.values
    if sub_cross.shape == (2, 2) and sub_cross.min() > 0:
        sub_t = Table2x2(sub_cross)
        labels.append(label)
        ors.append(sub_t.oddsratio)
        ci = sub_t.oddsratio_confint()
        lower_cis.append(ci[0])
        upper_cis.append(ci[1])

fig, ax = plt.subplots(figsize=(8, 5))
y_pos = range(len(labels))
ax.errorbar(ors, y_pos,
            xerr=[np.array(ors) - np.array(lower_cis), np.array(upper_cis) - np.array(ors)],
            fmt='D', color='black', capsize=3, markersize=5)
ax.axvline(x=1.0, color='gray', linestyle='--', linewidth=0.8)
ax.set_yticks(list(y_pos))
ax.set_yticklabels(labels)
ax.set_xlabel('Odds Ratio (95% CI)')
ax.set_xscale('log')
plt.tight_layout()
plt.savefig('forest_plot.png', dpi=150)
print(f'\nSaved forest_plot.png')

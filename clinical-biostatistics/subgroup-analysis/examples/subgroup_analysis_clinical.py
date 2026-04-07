'''Stratified subgroup analysis of a synthetic clinical trial dataset'''
# Reference: statsmodels 0.14+, scipy 1.12+, pandas 2.1+, matplotlib 3.8+ | Verify API if version differs

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.contingency_tables import StratifiedTable
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt

np.random.seed(42)
n = 600
df = pd.DataFrame({
    'treatment': np.random.choice(['Drug', 'Placebo'], n),
    'age_group': np.random.choice(['<65', '65+'], n, p=[0.6, 0.4]),
    'sex': np.random.choice(['M', 'F'], n),
    'site': np.random.choice(['Site_A', 'Site_B', 'Site_C'], n),
})

base_prob = 0.3
treatment_effect = -0.15
age_modifier = 0.1
df['prob'] = base_prob + (df['treatment'] == 'Drug').astype(float) * treatment_effect + (df['age_group'] == '65+').astype(float) * age_modifier
df['prob'] = df['prob'].clip(0.05, 0.95)
df['outcome'] = np.random.binomial(1, df['prob'])

tables = []
strata_names = []
for stratum in df['site'].unique():
    sub = df[df['site'] == stratum]
    t = pd.crosstab(sub['treatment'], sub['outcome']).values
    if t.shape == (2, 2):
        tables.append(t)
        strata_names.append(stratum)

st = StratifiedTable(tables)
print('=== Mantel-Haenszel Stratified Analysis (by site) ===')
print(st.summary())
print(f'\nPooled OR: {st.oddsratio_pooled:.3f}')
print(f'95% CI: {st.oddsratio_pooled_confint()}')

mh_test = st.test_null_odds()
print(f'MH test (H0: OR=1): statistic={mh_test.statistic:.3f}, p={mh_test.pvalue:.4f}')

bd_test = st.test_equal_odds()
print(f'Breslow-Day (H0: equal ORs): statistic={bd_test.statistic:.3f}, p={bd_test.pvalue:.4f}')

print('\n=== Interaction Test (treatment * age_group) ===')
interaction_model = smf.logit(
    'outcome ~ C(treatment, Treatment(reference="Placebo")) * C(age_group)', data=df
).fit(disp=0)
print(interaction_model.summary2().tables[1])

subgroups = ['age_group', 'sex', 'site']
labels, ors, lower_cis, upper_cis, pvals = [], [], [], [], []

for sg_var in subgroups:
    for group in sorted(df[sg_var].unique()):
        subset = df[df[sg_var] == group]
        sub_model = smf.logit(
            'outcome ~ C(treatment, Treatment(reference="Placebo"))', data=subset
        ).fit(disp=0)
        or_val = np.exp(sub_model.params.iloc[1])
        ci = np.exp(sub_model.conf_int().iloc[1])
        labels.append(f'{sg_var}: {group}')
        ors.append(or_val)
        lower_cis.append(ci[0])
        upper_cis.append(ci[1])
        pvals.append(sub_model.pvalues.iloc[1])

print('\n=== Subgroup-Specific ORs ===')
for i in range(len(labels)):
    print(f'{labels[i]:25s} OR={ors[i]:.3f} ({lower_cis[i]:.3f}-{upper_cis[i]:.3f}) p={pvals[i]:.4f}')

reject_holm, adjusted_holm, _, _ = multipletests(pvals, method='holm')
reject_fdr, adjusted_fdr, _, _ = multipletests(pvals, method='fdr_bh')

print('\n=== Multiplicity-Adjusted P-Values ===')
for i in range(len(labels)):
    print(f'{labels[i]:25s} Holm={adjusted_holm[i]:.4f}  BH={adjusted_fdr[i]:.4f}')

fig, ax = plt.subplots(figsize=(8, 6))
y_pos = range(len(labels))
ax.errorbar(ors, y_pos,
            xerr=[np.array(ors) - np.array(lower_cis), np.array(upper_cis) - np.array(ors)],
            fmt='D', color='black', capsize=3, markersize=5)
ax.axvline(x=1.0, color='gray', linestyle='--', linewidth=0.8)
ax.axvline(x=st.oddsratio_pooled, color='blue', linestyle=':', linewidth=0.8, alpha=0.5)
ax.set_yticks(y_pos)
ax.set_yticklabels(labels)
ax.set_xlabel('Odds Ratio (95% CI)')
ax.set_xscale('log')
ax.set_title('Subgroup Forest Plot')
plt.tight_layout()
plt.savefig('subgroup_forest_plot.png', dpi=150)
print('\nSaved subgroup_forest_plot.png')

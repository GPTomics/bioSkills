'''Chi-square and Fisher's exact tests for a clinical trial contingency table'''
# Reference: scipy 1.12+, statsmodels 0.14+, pandas 2.1+ | Verify API if version differs

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact
from statsmodels.stats.multitest import multipletests
from itertools import combinations

np.random.seed(42)
n = 200
df = pd.DataFrame({
    'treatment': np.random.choice(['Active', 'Placebo'], n),
    'outcome': np.random.choice(['Responder', 'Non-Responder', 'Partial'], n, p=[0.4, 0.35, 0.25])
})

table = pd.crosstab(df['treatment'], df['outcome'])
print('Contingency Table:')
print(table)

chi2, p, dof, expected = chi2_contingency(table, correction=False)
print(f'\nPearson Chi-Square: chi2={chi2:.4f}, dof={dof}, p={p:.4f}')

min_expected = expected.min()
print(f'Minimum expected count: {min_expected:.1f}')

if min_expected < 5:
    print('\nExpected count < 5 detected in RxC table, chi-square may be unreliable')
else:
    print('All expected counts >= 5, chi-square approximation is valid')

binary_table = pd.crosstab(
    df['treatment'],
    df['outcome'].map(lambda x: 'Responder' if x == 'Responder' else 'Non-Responder')
)
print(f'\n2x2 Table (Responder vs Non-Responder):')
print(binary_table)
chi2_2x2, p_2x2, _, expected_2x2 = chi2_contingency(binary_table, correction=False)
if expected_2x2.min() < 5:
    odds_ratio, p_fisher = fisher_exact(binary_table.values, alternative='two-sided')
    print(f'Fisher\'s exact: OR={odds_ratio:.4f}, p={p_fisher:.4f}')
else:
    print(f'Chi-square (2x2): chi2={chi2_2x2:.4f}, p={p_2x2:.4f}')

n_total = table.values.sum()
k = min(table.shape) - 1
cramers_v = np.sqrt(chi2 / (n_total * k))
print(f'\nCramer\'s V: {cramers_v:.4f}')

if k == 1:
    print(f'Interpretation: {"small" if cramers_v < 0.3 else "medium" if cramers_v < 0.5 else "large"} effect')
elif k == 2:
    print(f'Interpretation: {"small" if cramers_v < 0.21 else "medium" if cramers_v < 0.35 else "large"} effect')

if p < 0.05 and len(table.columns) > 2:
    print('\nPost-hoc pairwise comparisons (Holm correction):')
    categories = table.columns.tolist()
    pvalues = []
    comp_labels = []
    for cat1, cat2 in combinations(categories, 2):
        subset = df[df['outcome'].isin([cat1, cat2])]
        sub_table = pd.crosstab(subset['treatment'], subset['outcome'])
        _, p_val, _, _ = chi2_contingency(sub_table, correction=False)
        pvalues.append(p_val)
        comp_labels.append(f'{cat1} vs {cat2}')

    reject, adjusted_p, _, _ = multipletests(pvalues, method='holm')
    for label, raw_p, adj_p, sig in zip(comp_labels, pvalues, adjusted_p, reject):
        marker = '*' if sig else ''
        print(f'  {label}: raw p={raw_p:.4f}, adjusted p={adj_p:.4f} {marker}')
else:
    print('\nOverall test not significant or only 2 categories; post-hoc tests not needed')

results = pd.DataFrame({
    'statistic': ['chi2', 'dof', 'p_value', 'min_expected', 'cramers_v'],
    'value': [chi2, dof, p, min_expected, cramers_v]
})
results.to_csv('categorical_test_results.csv', index=False)
print('\nSaved categorical_test_results.csv')

'''Univariate differential metabolomics with explicit BH FDR and a volcano plot.

Demonstrates the two silent default traps: scipy ttest_ind defaults to Student
(equal_var=True) and statsmodels multipletests defaults to Holm-Sidak ('hs').
Both are overridden explicitly. Runs on synthetic data; the figure goes to a tempdir.
'''
# Reference: scipy 1.12+, statsmodels 0.14+, matplotlib 3.8+ | Verify API if version differs

import os
import tempfile
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

rng = np.random.default_rng(0)
n_per_group = 25
n_features = 400
n_true = 20                                   # only the first 20 features truly differ

samples = [f's{i}' for i in range(2 * n_per_group)]
features = [f'M{i}' for i in range(n_features)]
group = np.array(['control'] * n_per_group + ['case'] * n_per_group)

raw = rng.lognormal(mean=10, sigma=0.4, size=(n_features, 2 * n_per_group))
case_cols = np.where(group == 'case')[0]
raw[:n_true][:, case_cols] *= 2.0             # 2-fold up in case for the true features
intensities = pd.DataFrame(raw, index=features, columns=samples)

case = [s for s, g in zip(samples, group) if g == 'case']
ctrl = [s for s, g in zip(samples, group) if g == 'control']

logged = np.log2(intensities.replace(0, np.nan))   # transform before any parametric test

pvals, lfc = [], []
for feat in logged.index:
    a = logged.loc[feat, case].dropna().values
    b = logged.loc[feat, ctrl].dropna().values
    if len(a) >= 3 and len(b) >= 3:
        pvals.append(ttest_ind(a, b, equal_var=False)[1])   # Welch, not Student
        lfc.append(a.mean() - b.mean())                     # geometric-mean ratio on log scale
    else:
        pvals.append(np.nan)
        lfc.append(np.nan)

res = pd.DataFrame({'feature': logged.index, 'log2fc': lfc, 'pval': pvals}).dropna(subset=['pval'])
res['padj'] = multipletests(res['pval'], method='fdr_bh')[1]   # explicit BH, not 'hs'
res['hit'] = (res['padj'] < 0.05) & (res['log2fc'].abs() > 1)  # 2-fold + FDR 5%
res = res.sort_values('padj')

n_hit = int(res['hit'].sum())
true_recovered = res.loc[res['hit'], 'feature'].isin(features[:n_true]).sum()
print(f'Hits (padj<0.05 and |log2fc|>1): {n_hit} of {len(res)}')
print(f'True features recovered among hits: {true_recovered} of {n_true}')
print(res.head(10)[['feature', 'log2fc', 'pval', 'padj']].to_string(index=False))

plt.figure(figsize=(7, 5))
plt.scatter(res['log2fc'], -np.log10(res['pval']),
            c=np.where(res['hit'], 'firebrick', 'gray'), s=12, alpha=0.6)
plt.axhline(-np.log10(0.05), ls='--', color='black')
plt.axvline(1, ls='--', color='black')
plt.axvline(-1, ls='--', color='black')
plt.xlabel('log2 fold change')
plt.ylabel('-log10(p)')
out_path = os.path.join(tempfile.gettempdir(), 'metabolomics_volcano.png')
plt.savefig(out_path, dpi=120, bbox_inches='tight')
print(f'Volcano plot written to {out_path}')

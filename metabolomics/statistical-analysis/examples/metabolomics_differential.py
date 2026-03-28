'''Differential abundance analysis for metabolomics data with preprocessing'''
# Reference: scipy 1.12+, statsmodels 0.14+ | Verify API if version differs

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

intensities = pd.read_csv('feature_table.csv', index_col=0)
metadata = pd.read_csv('sample_info.csv')

groups = metadata.set_index('sample_id')['group']
group_levels = sorted(groups.unique())
ctrl_group, case_group = group_levels[0], group_levels[1]
ctrl_samples = groups[groups == ctrl_group].index.tolist()
case_samples = groups[groups == case_group].index.tolist()

# Preprocess: log2 transform (replace zeros with NaN to avoid -inf)
log2_data = np.log2(intensities.replace(0, np.nan))

# PQN normalization: reference spectrum -> quotients -> median quotient per sample
reference = log2_data.median(axis=1)
quotients = log2_data.div(reference, axis=0)
norm_factors = quotients.median(axis=0)
normalized = log2_data.div(norm_factors, axis=1)

# Welch's t-test per feature (equal_var=False is critical -- scipy defaults to Student's)
pvalues, log2fcs = [], []
for feature in normalized.index:
    case_vals = normalized.loc[feature, case_samples].dropna().values.astype(float)
    ctrl_vals = normalized.loc[feature, ctrl_samples].dropna().values.astype(float)
    if len(case_vals) >= 2 and len(ctrl_vals) >= 2:
        _, pval = ttest_ind(case_vals, ctrl_vals, equal_var=False)
        pvalues.append(pval)
        log2fcs.append(case_vals.mean() - ctrl_vals.mean())
    else:
        pvalues.append(np.nan)
        log2fcs.append(np.nan)

results = pd.DataFrame({'feature_id': normalized.index, 'log2fc': log2fcs, 'pvalue': pvalues})
results = results.dropna(subset=['pvalue'])

# BH FDR correction (statsmodels defaults to Holm-Sidak -- always specify fdr_bh)
_, results['padj'], _, _ = multipletests(results['pvalue'], method='fdr_bh')
results['significant'] = results['padj'] < 0.05

results = results.sort_values('padj')
results.to_csv('differential_results.csv', index=False)

n_sig = results['significant'].sum()
print(f'Significant features (padj < 0.05): {n_sig} of {len(results)}')
print(f'\nTop features by adjusted p-value:')
print(results.head(10)[['feature_id', 'log2fc', 'pvalue', 'padj']].to_string(index=False))

'''Differential abundance analysis for metabolomics data'''
# Reference: scipy 1.12+, statsmodels 0.14+ | Verify API if version differs

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

intensities = pd.read_csv('feature_table.csv', index_col=0)
metadata = pd.read_csv('sample_info.csv')

groups = metadata.set_index('sample_id')['group']
group_levels = sorted(groups.unique())
ctrl_group, treat_group = group_levels[0], group_levels[1]
ctrl_samples = groups[groups == ctrl_group].index.tolist()
treat_samples = groups[groups == treat_group].index.tolist()

# Log2 transform raw intensities; add pseudocount if zeros present
has_zeros = (intensities == 0).any().any()
log2_data = np.log2(intensities + 1) if has_zeros else np.log2(intensities)

# Welch's t-test per feature (equal_var=False is critical -- scipy defaults to Student's)
pvalues, log2fcs = [], []
for feature in log2_data.index:
    treat_vals = log2_data.loc[feature, treat_samples].values.astype(float)
    ctrl_vals = log2_data.loc[feature, ctrl_samples].values.astype(float)
    _, pval = ttest_ind(treat_vals, ctrl_vals, equal_var=False)
    pvalues.append(pval)
    log2fcs.append(treat_vals.mean() - ctrl_vals.mean())

results = pd.DataFrame({'feature_id': log2_data.index, 'log2fc': log2fcs, 'pvalue': pvalues})

# Benjamini-Hochberg FDR correction
_, padj, _, _ = multipletests(results['pvalue'], method='fdr_bh')
results['padj'] = padj
# FDR < 0.05 is the standard threshold for metabolomics differential abundance
results['significant'] = results['padj'] < 0.05

results = results.sort_values('padj')
results.to_csv('differential_results.csv', index=False)

n_sig = results['significant'].sum()
print(f'Significant features (padj < 0.05): {n_sig} of {len(results)}')
print(f'\nTop features by adjusted p-value:')
print(results.head(10)[['feature_id', 'log2fc', 'pvalue', 'padj']].to_string(index=False))

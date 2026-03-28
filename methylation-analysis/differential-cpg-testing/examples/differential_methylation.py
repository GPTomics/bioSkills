# Reference: scipy 1.12+, statsmodels 0.14+, pandas 2.2+, numpy 1.26+ | Verify API if version differs
# Per-CpG differential methylation: beta values, coverage filter, Welch's t-test, BH FDR

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

counts = pd.read_csv('bisulfite_counts.tsv', sep='\t', index_col='cpg_id')

case_samples = [c.replace('_meth', '') for c in counts.columns if c.startswith('case') and c.endswith('_meth')]
ctrl_samples = [c.replace('_meth', '') for c in counts.columns if c.startswith('ctrl') and c.endswith('_meth')]

case_betas = pd.DataFrame({s: counts[f'{s}_meth'] / counts[f'{s}_total'] for s in case_samples}, index=counts.index)
ctrl_betas = pd.DataFrame({s: counts[f'{s}_meth'] / counts[f'{s}_total'] for s in ctrl_samples}, index=counts.index)

# 10x minimum: provides 11 possible methylation levels per sample,
# adequate for reliable beta estimation and statistical power
MIN_COVERAGE = 10
total_cols = [c for c in counts.columns if c.endswith('_total')]
passes_coverage = (counts[total_cols] >= MIN_COVERAGE).all(axis=1)

case_betas = case_betas[passes_coverage]
ctrl_betas = ctrl_betas[passes_coverage]

# equal_var=False: Welch's t-test (scipy defaults to Student's with equal_var=True)
# Welch's is almost always more appropriate -- methylation variance differs between groups
t_stats, pvalues = ttest_ind(case_betas.values, ctrl_betas.values,
                              axis=1, equal_var=False, nan_policy='omit')

# method='fdr_bh': Benjamini-Hochberg FDR
# multipletests defaults to 'hs' (Holm-Sidak), not BH -- must specify explicitly
# 0.05: conventional FDR threshold for differential methylation
reject, padj, _, _ = multipletests(pvalues, alpha=0.05, method='fdr_bh')

mean_case = case_betas.mean(axis=1)
mean_ctrl = ctrl_betas.mean(axis=1)
delta_beta = mean_case - mean_ctrl

results = pd.DataFrame({
    'cpg_id': case_betas.index,
    'mean_case_beta': mean_case.values,
    'mean_ctrl_beta': mean_ctrl.values,
    'delta_beta': delta_beta.values,
    'pvalue': pvalues,
    'padj': padj,
    'significant': np.where(padj < 0.05, 'TRUE', 'FALSE')
})

results.to_csv('dmc.tsv', sep='\t', index=False)

n_sig = (results['significant'] == 'TRUE').sum()
print(f'CpGs tested: {len(results)}, significant (padj < 0.05): {n_sig}')

# Reference: numpy 1.26+, pandas 2.2+, scipy 1.12+, statsmodels 0.14+ | Verify API if version differs
# Differential protein abundance: log2 transform, median normalization, Welch's t-test, BH correction
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


def preprocess(intensities):
    log2_data = np.log2(intensities.replace(0, np.nan))
    sample_medians = log2_data.median(axis=0)
    global_median = sample_medians.median()
    return log2_data - sample_medians + global_median


def differential_abundance(normalized, case_cols, ctrl_cols):
    results = []
    for protein in normalized.index:
        case = normalized.loc[protein, case_cols].dropna()
        ctrl = normalized.loc[protein, ctrl_cols].dropna()
        if len(case) >= 2 and len(ctrl) >= 2:
            log2fc = case.mean() - ctrl.mean()
            _, pval = stats.ttest_ind(case, ctrl, equal_var=False)
            results.append({'protein_id': protein, 'log2fc': log2fc, 'pvalue': pval})

    df = pd.DataFrame(results)
    df['padj'] = multipletests(df['pvalue'], method='fdr_bh')[1]
    return df


intensities = pd.read_csv('protein_intensities.tsv', sep='\t', index_col=0)
groups = pd.read_csv('sample_groups.tsv', sep='\t')

case_cols = groups[groups['group'] == 'case']['sample_id'].tolist()
ctrl_cols = groups[groups['group'] == 'control']['sample_id'].tolist()

normalized = preprocess(intensities)
results = differential_abundance(normalized, case_cols, ctrl_cols)

results['significant'] = results['padj'].lt(0.05).map({True: 'TRUE', False: 'FALSE'})
results[['protein_id', 'log2fc', 'pvalue', 'padj', 'significant']].to_csv(
    'de_proteins.tsv', sep='\t', index=False)

n_sig = (results['significant'] == 'TRUE').sum()
print(f'Tested: {len(results)}, Significant: {n_sig}')

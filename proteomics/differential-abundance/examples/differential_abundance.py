# Reference: numpy 1.26+, pandas 2.2+, scipy 1.12+, statsmodels 0.14+ | Verify API if version differs
# Differential protein abundance: log2 transform, median normalization, Welch's t-test, BH correction.
# Welch + BH has NO variance moderation -- only appropriate at large n (>10/group). At n=3-5 use limma/DEqMS.
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

MIN_OBS_PER_GROUP = 2  # need >=2 non-missing values per group for a two-sample t with a variance estimate


def preprocess(intensities):
    log2_data = np.log2(intensities.replace(0, np.nan))  # zeros are undetected -> NaN to avoid -inf
    sample_medians = log2_data.median(axis=0)
    return log2_data - sample_medians + sample_medians.median()  # center every sample on the global median


def differential_abundance(normalized, case_cols, ctrl_cols):
    rows = []
    for protein in normalized.index:
        case, ctrl = normalized.loc[protein, case_cols].dropna(), normalized.loc[protein, ctrl_cols].dropna()
        if len(case) >= MIN_OBS_PER_GROUP and len(ctrl) >= MIN_OBS_PER_GROUP:
            _, pval = stats.ttest_ind(case, ctrl, equal_var=False)  # Welch; scipy defaults to Student's True
            rows.append({'protein': protein, 'log2fc': case.mean() - ctrl.mean(), 'pvalue': pval})
    df = pd.DataFrame(rows)
    df['padj'] = multipletests(df['pvalue'], method='fdr_bh')[1]  # default is Holm-Sidak; pass fdr_bh explicitly
    df['significant'] = df['padj'] < 0.05
    return df


def simulate_matrix(n_proteins=400, n_true=40, n_per_group=12, seed=0):
    rng = np.random.default_rng(seed)
    base = rng.normal(20, 0.4, size=(n_proteins, 2 * n_per_group))  # 0.4 log2 within-group spread (clean demo)
    base[:n_true, n_per_group:] += 1.5  # true up-regulation in case (1.5 log2 ~ 2.8-fold)
    cols = [f'ctrl_{i}' for i in range(n_per_group)] + [f'case_{i}' for i in range(n_per_group)]
    proteins = [f'P{i:04d}' for i in range(n_proteins)]
    return pd.DataFrame(2 ** base, index=proteins, columns=cols)


if __name__ == '__main__':
    intensities = simulate_matrix()
    ctrl_cols = [c for c in intensities.columns if c.startswith('ctrl')]
    case_cols = [c for c in intensities.columns if c.startswith('case')]

    normalized = preprocess(intensities)
    results = differential_abundance(normalized, case_cols, ctrl_cols)

    n_sig = int(results['significant'].sum())
    print(f'Tested: {len(results)}, Significant (padj<0.05): {n_sig}')
    print(results.sort_values('padj').head().to_string(index=False))

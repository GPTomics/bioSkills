# Reference: pandas 2.2+, numpy 1.26+, pydeseq2 0.5+ | Verify API if version differs

import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt


def log_cpm(counts, prior_count=2):
    '''Log2 CPM with scaled prior count (edgeR convention).'''
    lib_sizes = counts.sum(axis=0)
    cpm_vals = (counts + prior_count) / (lib_sizes + 2 * prior_count) * 1e6
    return np.log2(cpm_vals)


def counts_to_tpm(counts, gene_lengths):
    '''Convert raw counts to TPM. gene_lengths in bp.'''
    rate = counts.div(gene_lengths / 1000, axis=0)
    tpm = rate.div(rate.sum(axis=0), axis=1) * 1e6
    return tpm


def tmm_norm_factors(counts):
    '''Simplified TMM normalization factors.

    For production use, prefer edgeR::calcNormFactors via rpy2.
    '''
    lib_sizes = counts.sum(axis=0).values.astype(float)
    ref_idx = np.argmin(np.abs(np.log(lib_sizes) - np.mean(np.log(lib_sizes))))
    ref = counts.iloc[:, ref_idx].values.astype(float)
    factors = np.ones(counts.shape[1])

    for i in range(counts.shape[1]):
        if i == ref_idx:
            continue
        sample = counts.iloc[:, i].values.astype(float)
        keep = (ref > 0) & (sample > 0)
        if keep.sum() < 10:
            continue
        m = np.log2(sample[keep] / lib_sizes[i]) - np.log2(ref[keep] / lib_sizes[ref_idx])
        a = 0.5 * (np.log2(sample[keep] / lib_sizes[i]) + np.log2(ref[keep] / lib_sizes[ref_idx]))

        # Trim 30% of M-values and 5% of A-values
        m_lo, m_hi = np.percentile(m, [15, 85])
        a_lo, a_hi = np.percentile(a, [2.5, 97.5])
        keep2 = (m >= m_lo) & (m <= m_hi) & (a >= a_lo) & (a <= a_hi)
        if keep2.sum() > 0:
            factors[i] = 2 ** np.average(m[keep2])

    return factors / np.exp(np.mean(np.log(factors)))


def validate_size_factors(size_factors, sample_names):
    '''Check size factors for outliers.'''
    sf = pd.Series(size_factors, index=sample_names)
    print(f'Size factor range: {sf.min():.3f} - {sf.max():.3f}')
    print(f'Size factor median: {sf.median():.3f}')

    extreme = sf[(sf < 0.1) | (sf > 10)]
    if len(extreme) > 0:
        print(f'WARNING: {len(extreme)} samples have extreme size factors:')
        print(extreme)
        print('Investigate: contamination, degradation, or sample mislabeling')
    return sf


if __name__ == '__main__':
    np.random.seed(42)
    n_genes, n_samples = 5000, 6
    counts = pd.DataFrame(
        np.random.negative_binomial(5, 0.3, size=(n_genes, n_samples)),
        index=[f'gene_{i}' for i in range(n_genes)],
        columns=[f'sample_{i}' for i in range(n_samples)]
    )

    lcpm = log_cpm(counts, prior_count=2)
    print(f'log-CPM shape: {lcpm.shape}')
    print(f'log-CPM range: {lcpm.values.min():.2f} - {lcpm.values.max():.2f}')

    factors = tmm_norm_factors(counts)
    validate_size_factors(factors, counts.columns)

    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(lcpm.T)
    print(f'PC1 variance explained: {pca.explained_variance_ratio_[0]:.1%}')
    print(f'PC2 variance explained: {pca.explained_variance_ratio_[1]:.1%}')

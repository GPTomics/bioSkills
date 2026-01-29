'''Calculate translation efficiency from Ribo-seq and RNA-seq'''

import pandas as pd
import numpy as np

def load_counts(filepath):
    '''Load count matrix (genes x samples)'''
    return pd.read_csv(filepath, index_col=0)

def normalize_tpm(counts, lengths):
    '''Normalize counts to TPM

    TPM = (reads / gene_length_kb) / (sum_all_RPK / 1e6)
    Comparable across samples
    '''
    rpk = counts.div(lengths / 1000, axis=0)
    scaling = rpk.sum(axis=0) / 1e6
    tpm = rpk.div(scaling, axis=1)
    return tpm

def calculate_te(ribo_counts, rna_counts, gene_lengths):
    '''Calculate translation efficiency

    TE = Ribo-seq / RNA-seq (normalized)
    Higher TE = more efficiently translated

    Args:
        ribo_counts: Ribo-seq count matrix (genes x samples)
        rna_counts: RNA-seq count matrix (genes x samples)
        gene_lengths: Series with gene lengths

    Returns:
        DataFrame with log2 TE values
    '''
    # Normalize to TPM
    ribo_tpm = normalize_tpm(ribo_counts, gene_lengths)
    rna_tpm = normalize_tpm(rna_counts, gene_lengths)

    # Calculate TE with pseudocount
    # Pseudocount 0.1 TPM avoids log(0) and dampens noise in low-expressed genes
    PSEUDOCOUNT = 0.1
    te = (ribo_tpm + PSEUDOCOUNT) / (rna_tpm + PSEUDOCOUNT)

    # Log2 transform for symmetric fold changes
    log2_te = np.log2(te)

    return log2_te

def differential_te(te_matrix, condition_labels):
    '''Simple differential TE analysis

    Uses t-test for quick screening; riborex recommended for publication
    '''
    from scipy import stats

    conditions = pd.Series(condition_labels, index=te_matrix.columns)
    groups = conditions.unique()

    if len(groups) != 2:
        raise ValueError('Exactly 2 conditions required')

    group1_samples = conditions[conditions == groups[0]].index
    group2_samples = conditions[conditions == groups[1]].index

    results = []
    for gene in te_matrix.index:
        g1_values = te_matrix.loc[gene, group1_samples]
        g2_values = te_matrix.loc[gene, group2_samples]

        # Skip genes with no variance
        if g1_values.std() == 0 or g2_values.std() == 0:
            continue

        stat, pval = stats.ttest_ind(g1_values, g2_values)
        log2fc = g2_values.mean() - g1_values.mean()  # Already log2

        results.append({
            'gene': gene,
            'log2FC_TE': log2fc,
            'mean_group1': g1_values.mean(),
            'mean_group2': g2_values.mean(),
            'pvalue': pval
        })

    df = pd.DataFrame(results)

    # FDR correction using statsmodels for broader compatibility
    from statsmodels.stats.multitest import multipletests
    _, df['padj'], _, _ = multipletests(df['pvalue'], method='fdr_bh')

    return df.sort_values('padj')

# Example usage
if __name__ == '__main__':
    # Simulated example
    genes = [f'Gene{i}' for i in range(100)]
    samples = ['ctrl_1', 'ctrl_2', 'treat_1', 'treat_2']

    np.random.seed(42)

    # Simulate counts
    ribo_counts = pd.DataFrame(
        np.random.poisson(100, (100, 4)),
        index=genes, columns=samples
    )
    rna_counts = pd.DataFrame(
        np.random.poisson(500, (100, 4)),
        index=genes, columns=samples
    )
    gene_lengths = pd.Series(np.random.uniform(500, 5000, 100), index=genes)

    # Calculate TE
    log2_te = calculate_te(ribo_counts, rna_counts, gene_lengths)
    print('Translation Efficiency (log2):')
    print(log2_te.head())

    # Interpretation
    # log2 TE > 0: More ribosomes per mRNA than average
    # log2 TE < 0: Fewer ribosomes per mRNA than average
    # Delta log2 TE > 1: 2-fold increase in translation efficiency

    # Differential TE
    conditions = ['ctrl', 'ctrl', 'treat', 'treat']
    diff_te = differential_te(log2_te, conditions)
    print('\nDifferential TE (top genes):')
    print(diff_te.head(10))

    sig = diff_te[diff_te['padj'] < 0.05]
    print(f'\nGenes with significant TE change: {len(sig)}')

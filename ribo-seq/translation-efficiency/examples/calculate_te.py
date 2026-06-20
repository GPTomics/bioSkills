'''Per-gene translation efficiency (TE) as a ranking screen.

TE = RPF density / mRNA density, both counted over the CDS. This per-gene ratio
is for ranking and visualization ONLY. Differential TE is a condition-by-assay
interaction on raw counts and must use a count-based GLM (riborex / Xtail /
anota2seq / DESeq2 interaction) -- see the SKILL.md. anota2seq additionally
separates genuine translational control from buffering.
'''
# Reference: numpy 1.26+, pandas 2.2+, statsmodels 0.14+ | Verify API if version differs

import numpy as np
import pandas as pd

PSEUDOCOUNT = 0.1   # TPM pseudocount: avoids log(0), dampens low-count noise


def normalize_tpm(counts, cds_lengths):
    '''TPM over the CDS. cds_lengths in nucleotides, indexed like counts rows.'''
    rpk = counts.div(cds_lengths / 1000, axis=0)
    return rpk.div(rpk.sum(axis=0) / 1e6, axis=1)


def log2_te(ribo_counts, rna_counts, cds_lengths):
    '''Per-gene log2 TE matrix. Both assays counted over the CDS (not full transcript).'''
    ribo_tpm = normalize_tpm(ribo_counts, cds_lengths)
    rna_tpm = normalize_tpm(rna_counts, cds_lengths)
    return np.log2((ribo_tpm + PSEUDOCOUNT) / (rna_tpm + PSEUDOCOUNT))


def rank_te_change(te_matrix, conditions):
    '''Rank genes by mean log2 TE difference between two conditions (SCREEN ONLY).

    Reports a t-test for a quick look; it ignores count heteroskedasticity and is
    not a substitute for riborex/Xtail/anota2seq. Use only to prioritize.
    '''
    from scipy import stats
    from statsmodels.stats.multitest import multipletests

    cond = pd.Series(conditions, index=te_matrix.columns)
    groups = cond.unique()
    if len(groups) != 2:
        raise ValueError('Exactly 2 conditions required')
    g1 = cond[cond == groups[0]].index
    g2 = cond[cond == groups[1]].index

    rows = []
    for gene in te_matrix.index:
        a, b = te_matrix.loc[gene, g1], te_matrix.loc[gene, g2]
        if a.std() == 0 or b.std() == 0:
            continue
        _, p = stats.ttest_ind(a, b)
        rows.append({'gene': gene, 'log2FC_TE': b.mean() - a.mean(), 'pvalue': p})

    df = pd.DataFrame(rows)
    df['padj'] = multipletests(df['pvalue'], method='fdr_bh')[1]
    return df.sort_values('padj')


if __name__ == '__main__':
    rng = np.random.default_rng(42)
    genes = [f'Gene{i}' for i in range(100)]
    samples = ['ctrl_1', 'ctrl_2', 'treat_1', 'treat_2']
    ribo = pd.DataFrame(rng.poisson(100, (100, 4)), index=genes, columns=samples)
    rna = pd.DataFrame(rng.poisson(500, (100, 4)), index=genes, columns=samples)
    cds_len = pd.Series(rng.uniform(300, 5000, 100), index=genes)

    te = log2_te(ribo, rna, cds_len)
    print('Per-gene log2 TE (ranking screen only):')
    print(te.head())
    screen = rank_te_change(te, ['ctrl', 'ctrl', 'treat', 'treat'])
    print('\nTop TE-change candidates (verify with a count-based GLM):')
    print(screen.head())

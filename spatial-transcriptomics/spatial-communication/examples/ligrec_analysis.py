'''Ligand-receptor co-expression as a hypothesis, with honest multiple testing.

A ligrec hit is co-expression of ligand mRNA in a sender and receptor mRNA in a
receiver against a permutation null -- NOT signaling. squidpy ligrec is space-blind
(it permutes cluster labels); for distance-modeled inference use COMMOT (see SKILL.md).
'''
# Reference: squidpy 1.4+, scanpy 1.10+, statsmodels 0.14+ | Verify API if version differs

import squidpy as sq
from statsmodels.stats.multitest import multipletests

adata = sq.datasets.seqfish()
print(f'Loaded: {adata.n_obs} cells, {adata.obs["celltype_mapped_refined"].nunique()} cell types')

# CellPhoneDB permutation engine. threshold is the FRACTION of cells in a cluster that
# must express the gene (an expression floor), NOT a p-value cutoff. seqfish has no .raw.
res = sq.gr.ligrec(
    adata,
    cluster_key='celltype_mapped_refined',
    n_perms=100,
    threshold=0.01,
    use_raw=False,
    seed=0,
    n_jobs=1,
    show_progress_bar=False,
    copy=True,
)

pvalues = res['pvalues']
flat = pvalues.stack([0, 1], future_stack=True).rename('pval').reset_index().dropna(subset=['pval'])
flat['padj'] = multipletests(flat['pval'].values, method='fdr_bh')[1]
hits = flat[flat['padj'] < 0.05].sort_values('padj')

print(f'\n{len(hits)} co-expression hypotheses survive BH-FDR out of {len(flat)} tests')
print('These are candidates, not signaling: validate short-range hits against segmentation,')
print('then climb the confidence ladder (receiver-response -> protein co-localization -> perturbation).')
print('\nTop 10 by adjusted p-value:')
print(hits.head(10).to_string(index=False))

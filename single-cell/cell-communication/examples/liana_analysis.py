# Reference: liana-py 1.2+, matplotlib 3.8+, pandas 2.2+, scanpy 1.10+ | Verify API if version differs
import scanpy as sc
import liana as li
import pandas as pd
import matplotlib.pyplot as plt

adata = sc.read_h5ad('adata_annotated.h5ad')

li.mt.rank_aggregate(adata, groupby='cell_type', resource_name='consensus',
                     expr_prop=0.1, use_raw=False, n_perms=1000, verbose=True)

liana_results = adata.uns['liana_res']
print(f'Total interactions: {len(liana_results)}')

# rank_aggregate yields magnitude_rank and specificity_rank (no single 'liana_rank')
sig_interactions = liana_results[
    (liana_results['specificity_rank'] < 0.05) & (liana_results['magnitude_rank'] < 0.05)
].copy()
print(f'Robust interactions (both ranks < 0.05): {len(sig_interactions)}')

sig_interactions.to_csv('liana_significant_interactions.csv', index=False)

fig, ax = plt.subplots(figsize=(12, 8))
li.pl.dotplot(adata, colour='magnitude_rank', size='specificity_rank',
              top_n=50, orderby='magnitude_rank', orderby_ascending=True, ax=ax)
plt.tight_layout()
plt.savefig('liana_dotplot.png', dpi=150, bbox_inches='tight')
plt.close()

source_cells = ['Macrophage', 'Dendritic']
target_cells = ['T_cell', 'B_cell']
fig, ax = plt.subplots(figsize=(10, 8))
li.pl.dotplot(adata, colour='magnitude_rank', size='specificity_rank',
              source_labels=source_cells, target_labels=target_cells, ax=ax)
plt.tight_layout()
plt.savefig('liana_immune_crosstalk.png', dpi=150, bbox_inches='tight')
plt.close()

top_by_pair = sig_interactions.groupby(['source', 'target']).apply(
    lambda x: x.nsmallest(5, 'magnitude_rank')
).reset_index(drop=True)
top_by_pair.to_csv('liana_top_interactions_by_pair.csv', index=False)

adata.write('adata_with_liana.h5ad')

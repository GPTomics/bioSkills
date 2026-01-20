'''Ligand-receptor analysis for cell-cell communication'''

import squidpy as sq
import scanpy as sc
import pandas as pd

adata = sc.read_h5ad('clustered_spatial.h5ad')
print(f'Loaded: {adata.n_obs} cells/spots')
print(f'Cell types: {adata.obs["cell_type"].nunique()}')

sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=6)

print('\nRunning ligand-receptor analysis...')
sq.gr.ligrec(adata, cluster_key='cell_type', n_perms=100, threshold=0.01)

results = adata.uns['cell_type_ligrec']
pvalues = results['pvalues']
means = results['means']

sig_count = (pvalues < 0.05).sum().sum()
print(f'\nSignificant interactions (p < 0.05): {sig_count}')

interactions = []
for source_target in pvalues.index:
    for lr_pair in pvalues.columns:
        if pvalues.loc[source_target, lr_pair] < 0.05:
            source, target = source_target
            ligand, receptor = lr_pair
            interactions.append({
                'source': source, 'target': target,
                'ligand': ligand, 'receptor': receptor,
                'mean': means.loc[source_target, lr_pair],
                'pvalue': pvalues.loc[source_target, lr_pair],
            })

df = pd.DataFrame(interactions).sort_values('pvalue')
df.to_csv('significant_interactions.csv', index=False)
print(f'\nSaved {len(df)} interactions to significant_interactions.csv')
print('\nTop 10 interactions:')
print(df.head(10))

'''Compute spatial autocorrelation (Moran's I)'''

import squidpy as sq
import scanpy as sc

adata = sc.read_h5ad('preprocessed.h5ad')
print(f'Loaded: {adata.n_obs} spots')

# n_neighs=6: Default for hexagonal Visium spots; use 4 for square grids, adjust for spot density
sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=6)

# Limit to top 500 genes for speed; full genome takes longer but may reveal additional patterns
hvg = adata.var_names[adata.var['highly_variable']][:500]
print(f'\nComputing Moran\'s I for {len(hvg)} genes...')
sq.gr.spatial_autocorr(adata, mode='moran', genes=hvg)

results = adata.uns['moranI']
significant = results[results['pval_norm'] < 0.05].sort_values('I', ascending=False)

print(f'\nFound {len(significant)} spatially variable genes (p < 0.05)')
print('\nTop 10 spatially variable genes:')
print(significant.head(10)[['I', 'pval_norm']])

adata.write_h5ad('with_spatial_stats.h5ad')

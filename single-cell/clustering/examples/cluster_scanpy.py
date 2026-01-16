'''Cluster single-cell data with Scanpy'''

import scanpy as sc

adata = sc.read_h5ad('preprocessed.h5ad')
print(f'Loaded {adata.n_obs} cells')

sc.tl.pca(adata, n_comps=50)
sc.pl.pca_variance_ratio(adata, n_pcs=50, save='_variance.png')

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)

sc.tl.leiden(adata, resolution=0.5)
print(f'Found {adata.obs["leiden"].nunique()} clusters')
print(adata.obs['leiden'].value_counts())

sc.tl.umap(adata)
sc.pl.umap(adata, color='leiden', save='_clusters.png')

adata.write_h5ad('clustered.h5ad')
print('Saved clustered data')

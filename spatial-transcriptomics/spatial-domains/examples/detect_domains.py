'''Detect spatial domains with a BANKSY-style neighbor-augmented feature matrix.

Runs self-contained on the Squidpy Visium H&E dataset. Demonstrates the central
decisions: the expression-only salt-and-pepper baseline, the spatial-weight knob
(lambda) that controls over-smoothing, and a k+-1 sensitivity sweep. Heavy domain
packages (BayesSpace, STAGATE) are not required -- the neighbor-augmented matrix
reproduces the BANKSY idea with Squidpy + scanpy.
'''
# Reference: squidpy 1.4+, scanpy 1.10+, anndata 0.10+, scikit-learn 1.4+ | Verify API if version differs

import numpy as np
import scanpy as sc
import squidpy as sq
from sklearn.preprocessing import normalize
from sklearn.mixture import GaussianMixture

adata = sq.datasets.visium_hne_adata()
print(f'Loaded: {adata.n_obs} spots, {adata.n_vars} genes')

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=30)

# Expression-only baseline: spatially incoherent (salt-and-pepper).
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.leiden(adata, resolution=0.5, key_added='expr_leiden',
             flavor='igraph', n_iterations=2, directed=False)
print(f"Expression-only domains: {adata.obs['expr_leiden'].nunique()}")

# Spatial graph (Visium hex lattice -> 6 neighbors).
sq.gr.spatial_neighbors(adata, coord_type='grid', n_neighs=6)
W = normalize(adata.obsp['spatial_connectivities'], norm='l1', axis=1)

# Spatial-weight sweep: lambda is the over-smoothing knob (~0.8 for domains).
for lam in [0.2, 0.5, 0.8]:
    aug = np.concatenate(
        [np.sqrt(1 - lam) * adata.obsm['X_pca'],
         np.sqrt(lam) * (W @ adata.obsm['X_pca'])], axis=1)
    adata.obsm['X_banksy'] = aug
    sc.pp.neighbors(adata, use_rep='X_banksy', key_added='banksy')
    sc.tl.leiden(adata, resolution=0.5, key_added=f'domains_l{lam}',
                 neighbors_key='banksy', flavor='igraph', n_iterations=2, directed=False)
    print(f'lambda={lam}: {adata.obs[f"domains_l{lam}"].nunique()} domains')

# k as a biological choice -- report k+-1 instead of optimizing a score.
for k in [6, 7, 8]:
    gm = GaussianMixture(n_components=k, covariance_type='full', random_state=0)
    adata.obs[f'domains_k{k}'] = gm.fit_predict(adata.obsm['X_banksy']).astype(str)
print('k sweep done:', [f'k{k}' for k in [6, 7, 8]])

sc.tl.rank_genes_groups(adata, groupby='domains_l0.8', method='wilcoxon')
markers = sc.get.rank_genes_groups_df(adata, group=None)
print('\nTop markers per domain (lambda=0.8):')
print(markers.groupby('group').head(2)[['group', 'names', 'scores']])

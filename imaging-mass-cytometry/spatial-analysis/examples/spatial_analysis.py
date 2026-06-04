'''Per-image neighborhood enrichment with a pruned contact graph and a window-size sweep.

Encodes the two decisions the skill is built around: (1) run the permutation null
PER IMAGE and keep each z as a per-image summary (it is unbounded and graph-degree
dependent, so it is not comparable across images of different size), and (2) treat
cellular neighborhoods as exploratory by sweeping the window and checking stability.
'''
# Reference: squidpy 1.3+, scanpy 1.10+, anndata 0.10+, scikit-learn 1.4+, numpy 1.26+ | Verify API if version differs
import squidpy as sq
import anndata as ad
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

adata = ad.read_h5ad('imc_phenotyped.h5ad')
adata.obs['cell_type'] = adata.obs['cell_type'].astype('category')   # squidpy requires categorical

def pruned_contact_graph(sub, max_um=30.0):
    # Delaunay approximates physical adjacency, but invents long edges across lumen/necrosis;
    # prune anything longer than a biological contact distance
    sq.gr.spatial_neighbors(sub, coord_type='generic', delaunay=True)
    dist = sub.obsp['spatial_distances'].copy()
    dist.data[dist.data > max_um] = 0
    dist.eliminate_zeros()
    sub.obsp['spatial_connectivities'] = (dist > 0).astype(float)
    return sub

per_image_z = {}
for img_id, idx in adata.obs.groupby('image_id').groups.items():
    sub = pruned_contact_graph(adata[idx].copy())
    sq.gr.nhood_enrichment(sub, cluster_key='cell_type', seed=0)
    per_image_z[img_id] = sub.uns['cell_type_nhood_enrichment']['zscore']
print(f'Per-image z-scores for {len(per_image_z)} images; aggregate to PATIENT unit downstream')

def cellular_neighborhoods(adata, k_window, n_cn=10):
    sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=k_window)
    onehot = pd.get_dummies(adata.obs['cell_type']).values
    comp = adata.obsp['spatial_connectivities'] @ onehot
    comp = comp / comp.sum(axis=1, keepdims=True).clip(min=1)
    return KMeans(n_clusters=n_cn, random_state=0).fit_predict(np.asarray(comp))

# the window IS the spatial scale -- sweep it and report stability rather than trusting k=10
cn10 = cellular_neighborhoods(adata, k_window=10)
cn20 = cellular_neighborhoods(adata, k_window=20)
print(f'CN stability (ARI, window 10 vs 20): {adjusted_rand_score(cn10, cn20):.2f}  (low -> niches are scale-dependent)')
adata.obs['CN'] = cn20
adata.write('imc_spatial_analyzed.h5ad')

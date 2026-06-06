'''Why an out-of-distribution signal (not the softmax) catches novel cell types.

Runs end-to-end on synthetic latent coordinates: three reference cell types and
a query containing a fourth, NOVEL population that belongs to none of them. The
classifier softmax assigns every novel cell a confident label, while a distance-
based OOD signal (here, distance to the nearest reference cell -- the idea behind
Symphony's Mahalanobis gate) flags them. Weighted-kNN label-transfer uncertainty
(HLCA's 0.2 cutoff) is the other standard gate; it best catches cells that fall
BETWEEN reference types, while a distance gate best catches cells far from all of them.
'''
# Reference: scikit-learn 1.3+, anndata 0.10+ | Verify API if version differs

import numpy as np
import anndata as ad
from sklearn.neighbors import KNeighborsClassifier

rng = np.random.default_rng(0)

ref_centers = {'Tcell': [0, 0], 'Bcell': [6, 0], 'Myeloid': [0, 6]}
ref_latent = np.vstack([rng.normal(c, 0.6, (300, 2)) for c in ref_centers.values()])
ref_labels = np.repeat(list(ref_centers), 300)

query_known = np.vstack([rng.normal(c, 0.6, (80, 2)) for c in ref_centers.values()])
query_novel = rng.normal([12, 12], 0.6, (80, 2))      # hepatocyte-like: belongs to nothing
query_latent = np.vstack([query_known, query_novel])
novel_mask = np.array([False] * query_known.shape[0] + [True] * query_novel.shape[0])

knn = KNeighborsClassifier(n_neighbors=15).fit(ref_latent, ref_labels)
predicted = knn.predict(query_latent)
softmax_conf = knn.predict_proba(query_latent).max(axis=1)    # WRONG signal: which label, not whether

# Distance-based OOD gate: mean distance to the k nearest REFERENCE cells. Calibrate the
# threshold on the reference itself (99th percentile of its own nearest-neighbor distances).
ref_self_dist = knn.kneighbors(ref_latent)[0].mean(axis=1)
threshold = np.percentile(ref_self_dist, 99)
query_dist = knn.kneighbors(query_latent)[0].mean(axis=1)
gated = predicted.copy()
gated[query_dist > threshold] = 'Unknown'

print(f'Novel cells passing a 0.5 softmax filter (WRONG): {(softmax_conf[novel_mask] >= 0.5).mean():.0%}')
print(f'Novel cells caught by distance OOD gate (RIGHT):  {(gated[novel_mask] == "Unknown").mean():.0%}')
print(f'Known cells wrongly flagged Unknown:              {(gated[~novel_mask] == "Unknown").mean():.0%}')

adata_query = ad.AnnData(query_latent)
adata_query.obs['predicted_label'] = gated
adata_query.obs['ood_distance'] = query_dist
print(adata_query.obs['predicted_label'].value_counts())

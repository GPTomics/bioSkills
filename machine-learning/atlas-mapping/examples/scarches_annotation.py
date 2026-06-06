'''scArches scANVI label transfer with out-of-distribution gating (real-data template).

Requires a saved reference scANVI model and a query h5ad. The non-obvious point:
gate transferred labels on a weighted-kNN transfer uncertainty (does the cell
belong?), NOT on the classifier softmax max (which only says which label). See
ood_gating_demo.py for a synthetic, runnable version of the gating logic.
'''
# Reference: anndata 0.10+, scanpy 1.10+, scvi-tools 1.1+, scikit-learn 1.3+ | Verify API if version differs

import scvi
import scanpy as sc
import numpy as np
from sklearn.neighbors import KNeighborsClassifier

adata_ref = sc.read_h5ad('reference_labeled.h5ad')

# Reference scANVI head from a trained scVI model. unlabeled_category is REQUIRED.
scvi.model.SCVI.setup_anndata(adata_ref, layer='counts', batch_key='batch')
ref_vae = scvi.model.SCVI(adata_ref, n_latent=30, n_layers=2)
ref_vae.train(max_epochs=100, early_stopping=True)
ref_scanvi = scvi.model.SCANVI.from_scvi_model(ref_vae, unlabeled_category='Unknown', labels_key='cell_type')
ref_scanvi.train(max_epochs=20, n_samples_per_label=100)

adata_query = sc.read_h5ad('query.h5ad')
# Mandatory gene alignment: zero-pads missing, reorders. Silent corruption if skipped.
scvi.model.SCANVI.prepare_query_anndata(adata_query, ref_scanvi)
query_scanvi = scvi.model.SCANVI.load_query_data(adata_query, ref_scanvi)
# weight_decay=0.0 keeps the shared latent fixed so queries stay cross-comparable.
query_scanvi.train(max_epochs=100, plan_kwargs={'weight_decay': 0.0})

adata_query.obs['predicted_label'] = query_scanvi.predict()
adata_query.obsm['X_scANVI'] = query_scanvi.get_latent_representation()

# OOD gate on the shared latent. HLCA (Sikkema 2023) sets uncertainty > 0.2 to Unknown.
ref_latent = ref_scanvi.get_latent_representation()
knn = KNeighborsClassifier(n_neighbors=15, weights='distance').fit(ref_latent, adata_ref.obs['cell_type'])
uncertainty = 1.0 - knn.predict_proba(adata_query.obsm['X_scANVI']).max(axis=1)
adata_query.obs['transfer_uncertainty'] = uncertainty
adata_query.obs.loc[uncertainty > 0.2, 'predicted_label'] = 'Unknown'   # 0.2: HLCA default

print(f'Flagged Unknown (OOD / novel): {(uncertainty > 0.2).mean():.1%}')
print(adata_query.obs['predicted_label'].value_counts())

'''scArches transfer learning for single-cell annotation'''

import scvi
import scanpy as sc
import pandas as pd

adata_ref = sc.read_h5ad('reference.h5ad')
print(f'Reference: {adata_ref.n_obs} cells, {adata_ref.n_vars} genes')
print(f'Cell types: {adata_ref.obs["cell_type"].nunique()}')

# Train reference scVI model (skip if loading pre-trained)
scvi.model.SCVI.setup_anndata(adata_ref, layer='counts', batch_key='batch')
ref_vae = scvi.model.SCVI(adata_ref, n_latent=30, n_layers=2)
ref_vae.train(max_epochs=100, early_stopping=True)

# Convert to scANVI for label transfer
scvi.model.SCANVI.setup_anndata(adata_ref, layer='counts', batch_key='batch', labels_key='cell_type', unlabeled_category='Unknown')
ref_scanvi = scvi.model.SCANVI.from_scvi_model(ref_vae, labels_key='cell_type', unlabeled_category='Unknown')
ref_scanvi.train(max_epochs=50)
ref_scanvi.save('reference_scanvi/')
print('Saved reference scANVI model')

# Load query data
adata_query = sc.read_h5ad('query.h5ad')
print(f'\nQuery: {adata_query.n_obs} cells')

# Subset to reference genes (required)
common_genes = adata_ref.var_names.intersection(adata_query.var_names)
print(f'Common genes: {len(common_genes)}')
adata_query = adata_query[:, adata_ref.var_names].copy()

# Prepare and load query into scANVI
scvi.model.SCANVI.prepare_query_anndata(adata_query, ref_scanvi)
query_scanvi = scvi.model.SCANVI.load_query_data(adata_query, ref_scanvi)

# Surgical training: weight_decay=0.0 preserves reference model structure
# max_epochs=100: Usually sufficient for surgery
query_scanvi.train(max_epochs=100, plan_kwargs={'weight_decay': 0.0})

# Transfer labels
adata_query.obs['predicted_cell_type'] = query_scanvi.predict()
adata_query.obsm['X_scANVI'] = query_scanvi.get_latent_representation()

# Prediction confidence
soft_preds = query_scanvi.predict(soft=True)
adata_query.obs['prediction_confidence'] = soft_preds.max(axis=1)

# confidence < 0.5: Uncertain prediction, may be novel or poor quality
low_conf = adata_query.obs['prediction_confidence'] < 0.5
print(f'\nLow-confidence cells: {low_conf.sum()} ({low_conf.mean():.1%})')

print('\nPredicted cell type distribution:')
print(adata_query.obs['predicted_cell_type'].value_counts())

# Visualization
sc.pp.neighbors(adata_query, use_rep='X_scANVI')
sc.tl.umap(adata_query)
sc.pl.umap(adata_query, color=['predicted_cell_type', 'prediction_confidence'], save='_transferred_labels.png')

adata_query.write_h5ad('query_annotated.h5ad')
print('\nSaved annotated query to query_annotated.h5ad')

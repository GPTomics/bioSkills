'''Cell phenotyping from IMC data'''
import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd

# Load data
adata = ad.read_h5ad('imc_segmented.h5ad')
print(f'Loaded {adata.n_obs} cells, {adata.n_vars} markers')

# Arcsinh transform
adata.X = np.arcsinh(adata.X / 5)

# Phenotyping markers
markers = ['CD45', 'CD3', 'CD8', 'CD4', 'CD20', 'CD68', 'E-cadherin']
available = [m for m in markers if m in adata.var_names]
print(f'Using markers: {available}')

# Clustering
sc.pp.pca(adata, n_comps=min(15, len(available)))
sc.pp.neighbors(adata, n_neighbors=15)
sc.tl.leiden(adata, resolution=0.5)
sc.tl.umap(adata)

print(f'Found {adata.obs["leiden"].nunique()} clusters')

# Simple gating for major populations
def gate(marker, thresh):
    if marker in adata.var_names:
        return adata[:, marker].X.flatten() > thresh
    return np.zeros(adata.n_obs, dtype=bool)

adata.obs['CD45_pos'] = gate('CD45', 0.5)
adata.obs['CD3_pos'] = gate('CD3', 0.3)
adata.obs['CD8_pos'] = gate('CD8', 0.3)

# Assign basic types
def assign_type(row):
    if row['CD45_pos'] and row['CD3_pos'] and row['CD8_pos']:
        return 'CD8 T cell'
    elif row['CD45_pos'] and row['CD3_pos']:
        return 'T cell'
    elif row['CD45_pos']:
        return 'Immune (other)'
    return 'Non-immune'

adata.obs['cell_type'] = adata.obs.apply(assign_type, axis=1)

# Summary
print('\nCell type counts:')
print(adata.obs['cell_type'].value_counts())

# Save
adata.write('imc_phenotyped.h5ad')
print('\nSaved to imc_phenotyped.h5ad')

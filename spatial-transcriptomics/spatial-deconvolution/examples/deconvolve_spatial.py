'''Deconvolve spatial data with Tangram'''

import scanpy as sc
import tangram as tg
import numpy as np
import pandas as pd

adata_sc = sc.read_h5ad('reference_scrna.h5ad')
adata_sp = sc.read_h5ad('spatial_data.h5ad')

print(f'Reference: {adata_sc.n_obs} cells, {adata_sc.obs["cell_type"].nunique()} cell types')
print(f'Spatial: {adata_sp.n_obs} spots')

sc.pp.normalize_total(adata_sc)
sc.pp.log1p(adata_sc)

sc.tl.rank_genes_groups(adata_sc, groupby='cell_type', method='wilcoxon')
markers = sc.get.rank_genes_groups_df(adata_sc, group=None)
markers = markers[markers['pvals_adj'] < 0.01].groupby('group').head(100)
marker_genes = markers['names'].unique().tolist()
print(f'\nUsing {len(marker_genes)} marker genes')

tg.pp_adatas(adata_sc, adata_sp, genes=marker_genes)

print('\nMapping cells to space...')
ad_map = tg.map_cells_to_space(
    adata_sc, adata_sp,
    mode='clusters',
    cluster_label='cell_type',
    device='cpu',
)

tg.project_cell_annotations(ad_map, adata_sp, annotation='cell_type')

proportions = adata_sp.obsm['tangram_ct_pred']
cell_types = adata_sc.obs['cell_type'].cat.categories.tolist()

prop_df = pd.DataFrame(proportions, index=adata_sp.obs_names, columns=cell_types)
prop_df.to_csv('cell_type_proportions.csv')

adata_sp.obs['dominant_cell_type'] = cell_types[proportions.argmax(axis=1)]
adata_sp.write_h5ad('spatial_deconvolved.h5ad')

print('\nDeconvolution complete')
print(f'Saved: cell_type_proportions.csv, spatial_deconvolved.h5ad')

'''Preprocess single-cell data with Scanpy'''

import scanpy as sc

adata = sc.read_10x_mtx('filtered_feature_bc_matrix/', var_names='gene_symbols')
print(f'Raw: {adata.n_obs} cells, {adata.n_vars} genes')

adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs['pct_counts_mt'] < 20, :].copy()
print(f'Filtered: {adata.n_obs} cells, {adata.n_vars} genes')

adata.raw = adata.copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
print(f'HVGs: {adata.var.highly_variable.sum()}')

adata = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata, max_value=10)

adata.write_h5ad('preprocessed.h5ad')
print('Saved preprocessed data')

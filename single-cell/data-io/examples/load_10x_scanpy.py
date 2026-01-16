'''Load 10X Genomics data with Scanpy'''

import scanpy as sc

adata = sc.read_10x_mtx('filtered_feature_bc_matrix/', var_names='gene_symbols', cache=True)

print(f'Loaded {adata.n_obs} cells x {adata.n_vars} genes')
print(f'Cell IDs: {adata.obs_names[:5].tolist()}')
print(f'Gene names: {adata.var_names[:5].tolist()}')

adata.write_h5ad('pbmc.h5ad')
print('Saved to pbmc.h5ad')

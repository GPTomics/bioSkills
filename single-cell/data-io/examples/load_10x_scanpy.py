'''Load 10X Genomics data with Scanpy - multiple formats'''

import scanpy as sc

# Load from directory (MTX format - Cell Ranger output)
adata = sc.read_10x_mtx('filtered_feature_bc_matrix/', var_names='gene_symbols', cache=True)
print(f'Loaded {adata.n_obs} cells x {adata.n_vars} genes from MTX')

# Load from H5 file (Cell Ranger v3+ output)
adata_h5 = sc.read_10x_h5('filtered_feature_bc_matrix.h5')
print(f'Loaded {adata_h5.n_obs} cells x {adata_h5.n_vars} genes from H5')

# Load H5 with genome specification (for multi-species references)
adata_h5_genome = sc.read_10x_h5('filtered_feature_bc_matrix.h5', genome='GRCh38')

# Load H5 with gene IDs instead of symbols
adata_h5_ids = sc.read_10x_h5('filtered_feature_bc_matrix.h5', gex_only=True)

# Load raw counts (unfiltered) vs filtered
adata_raw = sc.read_10x_h5('raw_feature_bc_matrix.h5')
adata_filtered = sc.read_10x_h5('filtered_feature_bc_matrix.h5')
print(f'Raw: {adata_raw.n_obs} cells, Filtered: {adata_filtered.n_obs} cells')

# Save to h5ad format (Scanpy native)
adata.write_h5ad('pbmc.h5ad')
print('Saved to pbmc.h5ad')

'''Load 10X Genomics data with Scanpy - multiple formats'''
# Reference: cell ranger 8.0+, anndata 0.10+, numpy 1.26+, pandas 2.2+, scanpy 1.10+ | Verify API if version differs

import scanpy as sc

# Load from directory (MTX format - Cell Ranger output)
# cache=True keys the cache by PATH only: re-reading the same path after switching var_names returns the STALE object;
# delete the *.h5ad cache (or omit cache) when re-reading with different var_names
adata = sc.read_10x_mtx('filtered_feature_bc_matrix/', var_names='gene_symbols')
adata.var_names_make_unique()  # duplicate symbols get -1/-2 suffixes
print(f'Loaded {adata.n_obs} cells x {adata.n_vars} genes from MTX')

# Load from H5 file (Cell Ranger v3+ output)
adata_h5 = sc.read_10x_h5('filtered_feature_bc_matrix.h5')
print(f'Loaded {adata_h5.n_obs} cells x {adata_h5.n_vars} genes from H5')

# Load H5 with genome specification (for multi-species references)
adata_h5_genome = sc.read_10x_h5('filtered_feature_bc_matrix.h5', genome='GRCh38')

# Use stable Ensembl IDs instead of ambiguous symbols (MTX only; read_10x_h5 has no var_names arg, IDs land in var['gene_ids'])
adata_ids = sc.read_10x_mtx('filtered_feature_bc_matrix/', var_names='gene_ids')

# Keep non-GEX features (antibody capture, CRISPR guides); gex_only=True (default) drops them
adata_all = sc.read_10x_h5('filtered_feature_bc_matrix.h5', gex_only=False)

# Load raw counts (unfiltered) vs filtered
adata_raw = sc.read_10x_h5('raw_feature_bc_matrix.h5')
adata_filtered = sc.read_10x_h5('filtered_feature_bc_matrix.h5')
print(f'Raw: {adata_raw.n_obs} cells, Filtered: {adata_filtered.n_obs} cells')

# Save to h5ad format (Scanpy native)
adata.write_h5ad('pbmc.h5ad')
print('Saved to pbmc.h5ad')

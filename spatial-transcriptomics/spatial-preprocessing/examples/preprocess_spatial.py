'''Preprocess spatial transcriptomics data'''

import squidpy as sq
import scanpy as sc

adata = sq.read.visium('spaceranger_output/')
print(f'Loaded: {adata.n_obs} spots, {adata.n_vars} genes')

adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

print('\nQC metrics:')
print(f"  Total counts: {adata.obs['total_counts'].median():.0f} (median)")
print(f"  Genes/spot: {adata.obs['n_genes_by_counts'].median():.0f} (median)")
print(f"  MT%: {adata.obs['pct_counts_mt'].median():.1f}% (median)")

print(f'\nFiltering...')
sc.pp.filter_cells(adata, min_counts=1000)
sc.pp.filter_cells(adata, min_genes=500)
adata = adata[adata.obs['pct_counts_mt'] < 20].copy()
sc.pp.filter_genes(adata, min_cells=10)
print(f'After filtering: {adata.n_obs} spots, {adata.n_vars} genes')

adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat_v3', layer='counts')
print(f"HVGs: {adata.var['highly_variable'].sum()}")

sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50)

adata.write_h5ad('preprocessed.h5ad')
print('\nSaved to preprocessed.h5ad')

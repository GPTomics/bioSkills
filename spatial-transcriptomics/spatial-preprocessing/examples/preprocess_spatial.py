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
# Spatial data typically has higher background than scRNA-seq; thresholds depend on tissue
# min_counts=1000: Removes low-quality spots (tissue edges, empty areas). Adjust 500-2000 based on depth.
# min_genes=500: Ensures sufficient gene complexity. Lower for sparse tissues (muscle), higher for brain.
sc.pp.filter_cells(adata, min_counts=1000)
sc.pp.filter_cells(adata, min_genes=500)
# MT%<20: Standard cutoff. Spatial sections may have higher MT% due to tissue handling.
# For fragile tissues (e.g., adipose), consider MT%<30.
adata = adata[adata.obs['pct_counts_mt'] < 20].copy()
# min_cells=10: Gene must be in 10+ spots to be analyzed. Adjust for spot density.
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

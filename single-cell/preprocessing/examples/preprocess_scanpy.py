'''Preprocess single-cell data with Scanpy'''
# Reference: scanpy 1.10+ | Verify API if version differs

import scanpy as sc
import numpy as np
from scipy.stats import median_abs_deviation

adata = sc.read_10x_mtx('filtered_feature_bc_matrix/', var_names='gene_ids')  # gene_ids are stable across annotation releases
print(f'Raw: {adata.n_obs} cells, {adata.n_vars} genes')

adata.var['mt'] = adata.var['gene_symbols'].str.startswith('MT-')  # var_names are gene_ids here, so match symbols; mouse: 'mt-'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=[20], log1p=True, inplace=True)


def is_outlier(adata, metric, nmads):
    M = adata.obs[metric]
    return (M < np.median(M) - nmads * median_abs_deviation(M)) | (np.median(M) + nmads * median_abs_deviation(M) < M)


# MAD-adaptive thresholds, not fixed cutoffs: 5 MAD on counts/genes, 3 MAD on mito plus a hard cap
adata.obs['outlier'] = (is_outlier(adata, 'log1p_total_counts', 5) | is_outlier(adata, 'log1p_n_genes_by_counts', 5)
                        | is_outlier(adata, 'pct_counts_in_top_20_genes', 5))
# pct_counts_mt is a biology metric, not just quality: high baselines are normal in cardiomyocytes/hepatocytes/muscle
# and near-zero in nuclei. The 8% hard cap is tissue-dependent; raise or lower it for the tissue being analyzed
adata.obs['mt_outlier'] = is_outlier(adata, 'pct_counts_mt', 3) | (adata.obs['pct_counts_mt'] > 8)
adata = adata[~(adata.obs['outlier'] | adata.obs['mt_outlier'])].copy()
sc.pp.filter_genes(adata, min_cells=3)  # drop genes detected in fewer than 3 cells
print(f'Filtered: {adata.n_obs} cells, {adata.n_vars} genes')

adata.layers['counts'] = adata.X.copy()  # stash raw counts: HVG (seurat_v3) and downstream count models need them
adata.raw = adata

sc.pp.normalize_total(adata)  # target_sum=None scales each cell to the dataset median; 1e4 ('CP10k') is an arbitrary legacy choice
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat_v3', layer='counts')  # seurat_v3 reads RAW counts
print(f'HVGs: {adata.var.highly_variable.sum()}')

sc.pp.pca(adata, mask_var='highly_variable')  # PCA on log-normalized HVG; scaling is increasingly skipped

adata.write_h5ad('preprocessed.h5ad')
print('Saved preprocessed data')

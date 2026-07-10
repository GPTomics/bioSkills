# Reference: matplotlib 3.8+, numpy 1.26+, scanpy 1.10+, squidpy 1.3+ | Verify API if version differs
# Complete Visium spatial transcriptomics workflow

import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
import numpy as np
import os

sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=100, facecolor='white')

# Configuration
data_dir = 'spaceranger_output'
output_dir = 'visium_results'
os.makedirs(output_dir, exist_ok=True)
os.makedirs(f'{output_dir}/plots', exist_ok=True)

# === Step 1: Load Data ===
print('=== Step 1: Loading Data ===')
adata = sq.read.visium(data_dir)
adata.var_names_make_unique()
print(f'Loaded: {adata.n_obs} spots, {adata.n_vars} genes')

# === Step 2: QC ===
print('=== Step 2: Quality Control ===')
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], inplace=True)

# QC plots
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
sc.pl.spatial(adata, color='total_counts', ax=axes[0, 0], show=False)
sc.pl.spatial(adata, color='n_genes_by_counts', ax=axes[0, 1], show=False)
sc.pl.spatial(adata, color='pct_counts_mt', ax=axes[0, 2], show=False)
# sc.pl.violin takes a single Axes, not an array -> one violin per axis
for a, key in zip(axes[1, :], ['total_counts', 'n_genes_by_counts', 'pct_counts_mt']):
    sc.pl.violin(adata, key, jitter=0.4, ax=a, show=False)
plt.tight_layout()
plt.savefig(f'{output_dir}/plots/qc_metrics.pdf')
plt.close()

# Filter. These floors are VISIUM defaults (a spot is a 1-10-cell mixture). For
# IMAGING data (Xenium/MERFISH) use ~10 transcripts/cell, NOT 500, or this deletes
# nearly every real cell; mito QC is usually impossible (mito off-panel).
print(f'Before filtering: {adata.n_obs} spots')
sc.pp.filter_cells(adata, min_counts=500)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=10)
adata = adata[adata.obs.pct_counts_mt < 25, :]
print(f'After filtering: {adata.n_obs} spots')

# === Step 3: Normalization ===
print('=== Step 3: Normalization ===')
adata.layers['counts'] = adata.X.copy()
# Library size partly carries biology in spatial data (cells-per-spot); total-count
# normalization is a Visium starting point, not a universal default (see spatial-preprocessing).
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
print(f'HVGs: {adata.var.highly_variable.sum()}')

# === Step 4: Dimensionality Reduction & Clustering ===
print('=== Step 4: Clustering ===')
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5, flavor='igraph', n_iterations=2, directed=False)

# Cluster visualization
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
sc.pl.umap(adata, color='leiden', ax=axes[0], show=False)
sc.pl.spatial(adata, color='leiden', spot_size=1.5, ax=axes[1], show=False)
plt.tight_layout()
plt.savefig(f'{output_dir}/plots/clusters.pdf')
plt.close()

print(f'Clusters: {adata.obs["leiden"].nunique()}')

# === Step 5: Spatial Analysis ===
print('=== Step 5: Spatial Analysis ===')

# Spatial neighbors. Visium is a hex lattice -> coord_type='grid' (n_neighs=6);
# 'generic' kNN is for imaging point clouds (see spatial-neighbors).
sq.gr.spatial_neighbors(adata, coord_type='grid', n_neighs=6)

# Neighborhood enrichment
sq.gr.nhood_enrichment(adata, cluster_key='leiden')
sq.pl.nhood_enrichment(adata, cluster_key='leiden')
plt.savefig(f'{output_dir}/plots/nhood_enrichment.pdf')
plt.close()

# Co-occurrence
sq.gr.co_occurrence(adata, cluster_key='leiden')
sq.pl.co_occurrence(adata, cluster_key='leiden', clusters=['0', '1'])
plt.savefig(f'{output_dir}/plots/co_occurrence.pdf')
plt.close()

# Spatially variable genes. Gate on FDR (not raw I); and a top-Moran gene is usually
# a marker of a spatially-clustered cell TYPE (composition), not a gene regulated
# WITHIN a type -- intersect with non-HVG for the latter (see spatial-statistics).
print('Finding spatially variable genes...')
sq.gr.spatial_autocorr(adata, mode='moran', n_perms=100, n_jobs=4)
moran = adata.uns['moranI']
svg = moran[moran['pval_norm_fdr_bh'] < 0.05].sort_values('I', ascending=False)
top_svg = svg.head(20).index.tolist()
print(f'Top SVGs (FDR<0.05): {top_svg[:5]}')

# Plot top SVGs
sc.pl.spatial(adata, color=top_svg[:4], ncols=2, spot_size=1.5, cmap='viridis')
plt.savefig(f'{output_dir}/plots/top_svg.pdf')
plt.close()

# === Step 6: Cluster Markers ===
# On Visium these are markers of spot REGIONS/niches (mixtures), not pure cell types;
# for cell-type composition deconvolve first (see spatial-deconvolution).
print('=== Step 6: Cluster Markers ===')
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
markers = sc.get.rank_genes_groups_df(adata, group=None)
markers.to_csv(f'{output_dir}/cluster_markers.csv', index=False)

# Marker plots
# scanpy save= appends onto ./figures/dotplot_ ; use a plain suffix (not a path)
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, save='_markers_dotplot.pdf')

# === Step 7: Save Results ===
print('=== Step 7: Saving Results ===')
svg.to_csv(f'{output_dir}/spatially_variable_genes.csv')
adata.write(f'{output_dir}/visium_analyzed.h5ad')

print(f'\n=== Analysis Complete ===')
print(f'Results saved to: {output_dir}/')
print(f'  - Processed data: visium_analyzed.h5ad')
print(f'  - SVGs: spatially_variable_genes.csv')
print(f'  - Markers: cluster_markers.csv')
print(f'  - Plots: plots/')

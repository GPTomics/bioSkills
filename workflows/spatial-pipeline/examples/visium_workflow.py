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
sc.pl.violin(adata, ['total_counts', 'n_genes_by_counts', 'pct_counts_mt'],
             jitter=0.4, ax=axes[1, :], show=False)
plt.tight_layout()
plt.savefig(f'{output_dir}/plots/qc_metrics.pdf')
plt.close()

# Filter
print(f'Before filtering: {adata.n_obs} spots')
sc.pp.filter_cells(adata, min_counts=500)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=10)
adata = adata[adata.obs.pct_counts_mt < 25, :]
print(f'After filtering: {adata.n_obs} spots')

# === Step 3: Normalization ===
print('=== Step 3: Normalization ===')
adata.layers['counts'] = adata.X.copy()
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
sc.tl.leiden(adata, resolution=0.5)

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

# Spatial neighbors
sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=6)

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

# Spatially variable genes
print('Finding spatially variable genes...')
sq.gr.spatial_autocorr(adata, mode='moran', n_perms=100, n_jobs=4)
svg = adata.uns['moranI'].sort_values('I', ascending=False)
top_svg = svg.head(20).index.tolist()
print(f'Top SVGs: {top_svg[:5]}')

# Plot top SVGs
sc.pl.spatial(adata, color=top_svg[:4], ncols=2, spot_size=1.5, cmap='viridis')
plt.savefig(f'{output_dir}/plots/top_svg.pdf')
plt.close()

# === Step 6: Cluster Markers ===
print('=== Step 6: Cluster Markers ===')
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
markers = sc.get.rank_genes_groups_df(adata, group=None)
markers.to_csv(f'{output_dir}/cluster_markers.csv', index=False)

# Marker plots
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, save=f'_{output_dir}/plots/markers_dotplot.pdf')

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

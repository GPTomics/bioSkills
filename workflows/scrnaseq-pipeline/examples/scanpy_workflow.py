# Reference: matplotlib 3.8+, numpy 1.26+, pandas 2.2+, scanpy 1.10+ | Verify API if version differs
# Complete single-cell RNA-seq workflow with Scanpy

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=100, facecolor='white')

# Public scRNA-seq datasets:
# - 10x Genomics: https://www.10xgenomics.com/datasets (PBMC 3k, 10k datasets)
# - CELLxGENE: https://cellxgene.cziscience.com (curated annotated datasets)
# - GEO: GSE149173 (COVID-19 PBMC), GSE136831 (human lung)
# - Human Cell Atlas: https://data.humancellatlas.org
# - Scanpy built-in: sc.datasets.pbmc3k_processed()

# Configuration
data_path = 'filtered_feature_bc_matrix.h5'
output_dir = 'scrnaseq_results_scanpy'
os.makedirs(output_dir, exist_ok=True)
os.makedirs(f'{output_dir}/plots', exist_ok=True)

# === Step 1: Load Data ===
print('Loading data...')
adata = sc.read_10x_h5(data_path)
adata.var_names_make_unique()
print(f'Initial: {adata.n_obs} cells, {adata.n_vars} genes')

# === Step 2: QC Metrics ===
print('Calculating QC metrics...')
adata.var['mt'] = adata.var_names.str.startswith('MT-')
adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)

# QC plots (one violin per axis; sc.pl.violin takes a single Axes, not an array)
fig, axes = plt.subplots(1, 4, figsize=(16, 4))
for a, key in zip(axes[:3], ['n_genes_by_counts', 'total_counts', 'pct_counts_mt']):
    sc.pl.violin(adata, key, jitter=0.4, ax=a, show=False)
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', color='pct_counts_mt', ax=axes[3], show=False)
plt.savefig(f'{output_dir}/plots/qc_metrics.pdf')
plt.close()

# Filter cells
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs.n_genes_by_counts < 5000, :]
adata = adata[adata.obs.pct_counts_mt < 20, :]
print(f'After QC: {adata.n_obs} cells')

# === Step 3: Doublet Detection ===
print('Detecting doublets...')
# Expected rate from recovered cells (~0.8% per 1,000), not the 0.05 placeholder default
expected_rate = 0.008 * adata.n_obs / 1000
sc.pp.scrublet(adata, expected_doublet_rate=expected_rate)
doublet_rate = adata.obs['predicted_doublet'].sum() / len(adata)
print(f'Doublet rate: {doublet_rate:.1%}')

adata = adata[~adata.obs['predicted_doublet'], :]
print(f'After doublet removal: {adata.n_obs} cells')

# === Step 4: Normalization ===
print('Normalizing...')
adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# HVGs
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
print(f'Highly variable genes: {adata.var.highly_variable.sum()}')

# sc.pl.highly_variable_genes has no ax= param; let scanpy manage the figure, then save.
sc.pl.highly_variable_genes(adata, show=False)
plt.savefig(f'{output_dir}/plots/hvgs.pdf')
plt.close()

# === Step 5: Dimensionality Reduction ===
print('Running PCA...')
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50)

# Elbow plot (pca_variance_ratio has no ax= param)
sc.pl.pca_variance_ratio(adata, n_pcs=50, show=False)
plt.savefig(f'{output_dir}/plots/elbow.pdf')
plt.close()

print('Computing neighbors and UMAP...')
n_pcs = 30
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=n_pcs)
sc.tl.umap(adata)

# === Step 6: Clustering ===
print('Clustering...')
# Pin flavor='igraph' for reproducibility (scanpy 1.10 default flips to igraph)
for res in [0.2, 0.4, 0.5, 0.6, 0.8, 1.0]:
    sc.tl.leiden(adata, resolution=res, key_added=f'leiden_{res}',
                 flavor='igraph', n_iterations=2, directed=False)

# UMAP plots
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
for ax, res in zip(axes.flat, [0.2, 0.4, 0.6, 0.8]):
    sc.pl.umap(adata, color=f'leiden_{res}', ax=ax, show=False, title=f'res={res}')
plt.tight_layout()
plt.savefig(f'{output_dir}/plots/umap_resolutions.pdf')
plt.close()

# Set default clustering
adata.obs['leiden'] = adata.obs['leiden_0.5']
print(f'Clusters (res=0.5): {adata.obs["leiden"].nunique()}')

# === Step 7: Find Markers ===
print('Finding marker genes...')
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

# Plot markers (rank_genes_groups draws its own per-group panel grid; a pre-made ax would sit empty behind it)
sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False, show=False)
plt.savefig(f'{output_dir}/plots/markers.pdf')
plt.close()

# Export markers
markers = sc.get.rank_genes_groups_df(adata, group=None)
markers.to_csv(f'{output_dir}/all_markers.csv', index=False)

# Top markers per cluster
top_markers = markers.groupby('group').head(10)
top_markers.to_csv(f'{output_dir}/top10_markers.csv', index=False)

# Heatmap: sc.pl.heatmap builds its own multi-axis gridspec, so it takes NO ax= (unlike violin);
# pass show=False + save= and let scanpy write it.
top_genes = markers.groupby('group').head(5)['names'].tolist()
sc.pl.heatmap(adata, var_names=top_genes[:50], groupby='leiden', show=False,
              save='_marker_heatmap.pdf')   # writes figures/heatmap_marker_heatmap.pdf

# === Step 8: Save Results ===
print('Saving results...')
adata.write(f'{output_dir}/adata.h5ad')

# Summary
print('\n=== Pipeline Complete ===')
print(f'Final cells: {adata.n_obs}')
print(f'Clusters: {adata.obs["leiden"].nunique()}')
print(f'Results saved to: {output_dir}')
print('\nTop 3 markers per cluster:')
print(markers.groupby('group').head(3)[['group', 'names', 'logfoldchanges', 'pvals_adj']])

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

# QC plots
fig, axes = plt.subplots(1, 4, figsize=(16, 4))
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], jitter=0.4, ax=axes[:3], show=False)
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
sc.pp.scrublet(adata)
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

fig, ax = plt.subplots(figsize=(8, 6))
sc.pl.highly_variable_genes(adata, ax=ax, show=False)
plt.savefig(f'{output_dir}/plots/hvgs.pdf')
plt.close()

# === Step 5: Dimensionality Reduction ===
print('Running PCA...')
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50)

# Elbow plot
fig, ax = plt.subplots(figsize=(6, 4))
sc.pl.pca_variance_ratio(adata, n_pcs=50, ax=ax, show=False)
plt.savefig(f'{output_dir}/plots/elbow.pdf')
plt.close()

print('Computing neighbors and UMAP...')
n_pcs = 30
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=n_pcs)
sc.tl.umap(adata)

# === Step 6: Clustering ===
print('Clustering...')
for res in [0.2, 0.4, 0.6, 0.8, 1.0]:
    sc.tl.leiden(adata, resolution=res, key_added=f'leiden_{res}')

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

# Plot markers
fig, ax = plt.subplots(figsize=(12, 8))
sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False, ax=ax, show=False)
plt.savefig(f'{output_dir}/plots/markers.pdf')
plt.close()

# Export markers
markers = sc.get.rank_genes_groups_df(adata, group=None)
markers.to_csv(f'{output_dir}/all_markers.csv', index=False)

# Top markers per cluster
top_markers = markers.groupby('group').head(10)
top_markers.to_csv(f'{output_dir}/top10_markers.csv', index=False)

# Heatmap
top_genes = markers.groupby('group').head(5)['names'].tolist()
fig, ax = plt.subplots(figsize=(12, 10))
sc.pl.heatmap(adata, var_names=top_genes[:50], groupby='leiden', ax=ax, show=False)
plt.savefig(f'{output_dir}/plots/marker_heatmap.pdf')
plt.close()

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

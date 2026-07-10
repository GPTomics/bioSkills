"""Reference: scanpy 1.10+, umap-learn 0.5+, openTSNE 1.0+, phate 1.0+ | Verify API if version differs

PhD-level dimensionality-reduction with explicit hyperparameters, reproducible seeds, and
Chari-Pachter caveat communicated. Encodes four traps: random seed, PCA initialization for
t-SNE (Kobak-Berens 2019), scanpy save path trap, variance-labeled PCA axes.
"""
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import umap
import openTSNE
import phate

# Default DPI for publication; figdir explicit
sc.set_figure_params(dpi_save=300, figsize=(4, 4), frameon=False)
sc.settings.figdir = './figures/'                # REPLACE WITH ABSOLUTE PATH IN PRODUCTION; relative ./figures/ for testing

adata = sc.read_h5ad('processed.h5ad')           # QC-filtered counts; the category README covers upstream QC

# 1. PCA -- variance-explained axes, deterministic
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor='seurat_v3')
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50, random_state=42)

var = adata.uns['pca']['variance_ratio']
fig, ax = plt.subplots(figsize=(4, 4))
ax.scatter(adata.obsm['X_pca'][:, 0], adata.obsm['X_pca'][:, 1],
           c=adata.obs['condition'].astype('category').cat.codes,
           cmap='tab10', alpha=0.7, s=10, rasterized=True)
ax.set_xlabel(f'PC1 ({var[0]*100:.1f}%)')                # variance-labeled
ax.set_ylabel(f'PC2 ({var[1]*100:.1f}%)')
ax.set_title('PCA')
plt.savefig('pca.pdf', bbox_inches='tight', dpi=300)

# Scree plot for choosing n_pcs downstream
fig, ax = plt.subplots(figsize=(5, 3))
ax.plot(range(1, 51), var, 'o-', markersize=3)
ax.set_xlabel('Principal component'); ax.set_ylabel('Variance explained')

# 2. UMAP -- standard scanpy workflow, explicit hyperparameters
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=50)
sc.tl.umap(adata, min_dist=0.3, random_state=42)
sc.tl.leiden(adata, resolution=0.5, random_state=42)
sc.pl.umap(adata, color='leiden', palette='tab20',
           legend_loc='on data', legend_fontsize=7,
           save='_clusters.pdf')                          # writes to figdir/umap_clusters.pdf

# 3. t-SNE -- Kobak-Berens defaults (PCA init, learning_rate=n/12, perplexity 30)
X_pca = adata.obsm['X_pca'][:, :50]
n = X_pca.shape[0]
tsne_embedding = openTSNE.TSNE(
    perplexity=30,
    n_iter=750,
    initialization='pca',                                  # critical -- not 'random'
    learning_rate=n / 12,                                  # Kobak-Berens
    n_jobs=-1,
    random_state=42).fit(X_pca)

fig, ax = plt.subplots(figsize=(4, 4))
ax.scatter(tsne_embedding[:, 0], tsne_embedding[:, 1],
           c=adata.obs['leiden'].astype('category').cat.codes,
           cmap='tab20', s=2, alpha=0.7, rasterized=True)
ax.set_xlabel('t-SNE 1'); ax.set_ylabel('t-SNE 2')         # no variance units
ax.set_title(f't-SNE (perp=30, lr=n/12, PCA init)')
plt.savefig('tsne.pdf', bbox_inches='tight', dpi=300)

# 4. PHATE -- for continuous trajectories (not for discrete clusters)
phate_op = phate.PHATE(knn=10, decay=40, t='auto', n_jobs=-1, random_state=42)
phate_emb = phate_op.fit_transform(X_pca)

fig, ax = plt.subplots(figsize=(4, 4))
ax.scatter(phate_emb[:, 0], phate_emb[:, 1],
           c=adata.obs['pseudotime'] if 'pseudotime' in adata.obs else None,
           cmap='viridis', s=3, alpha=0.7, rasterized=True)
ax.set_xlabel('PHATE 1'); ax.set_ylabel('PHATE 2')
ax.set_title('PHATE')

# 5. CAPTION TEMPLATE for any 2D embedding figure
caption = ('UMAP (n_neighbors=30, min_dist=0.3, random_state=42) of N={n} cells. '
           'Colored by Leiden cluster (resolution=0.5). '
           '2D embeddings preserve local neighborhoods only; distances between clusters '
           'and density within clusters are NOT biological signals (Chari & Pachter 2023).')

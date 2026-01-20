import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
import squidpy as sq
from pathlib import Path

# === CONFIGURATION ===
data_dir = Path('steinbock_output')
output_dir = Path('results')
output_dir.mkdir(exist_ok=True)

# === 1. LOAD DATA ===
print('Loading steinbock outputs...')
intensities = pd.read_csv(data_dir / 'intensities.csv', index_col=0)
regionprops = pd.read_csv(data_dir / 'regionprops.csv', index_col=0)

print(f'Loaded {len(intensities)} cells, {len(intensities.columns)} markers')

# === 2. CREATE ANNDATA ===
adata = ad.AnnData(X=intensities.values, obs=regionprops.copy(),
                   var=pd.DataFrame(index=intensities.columns))
adata.obs['cell_id'] = intensities.index.values
adata.obs['image_id'] = [idx.rsplit('_', 1)[0] for idx in intensities.index]
adata.obsm['spatial'] = regionprops[['centroid_y', 'centroid_x']].values

# === 3. PREPROCESSING ===
print('Preprocessing...')
adata.X = np.arcsinh(adata.X / 5)
adata.raw = adata.copy()
sc.pp.scale(adata, max_value=10)

# === 4. DIMENSIONALITY REDUCTION ===
print('Running dimensionality reduction...')
sc.pp.pca(adata, n_comps=20)
sc.pp.neighbors(adata, n_neighbors=15)
sc.tl.umap(adata)

# === 5. CLUSTERING ===
print('Clustering...')
sc.tl.leiden(adata, resolution=0.8)
print(f'Found {adata.obs["leiden"].nunique()} clusters')

# === 6. MARKER ANALYSIS ===
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

# === 7. SPATIAL ANALYSIS ===
print('Running spatial analysis...')
sq.gr.spatial_neighbors(adata, coord_type='generic', delaunay=True)
sq.gr.nhood_enrichment(adata, cluster_key='leiden')

# === 8. VISUALIZATION ===
print('Generating plots...')
sc.pl.umap(adata, color=['leiden'], save='_clusters.png')
sq.pl.nhood_enrichment(adata, cluster_key='leiden')

# Spatial plot for first image
first_image = adata.obs['image_id'].iloc[0]
adata_img = adata[adata.obs['image_id'] == first_image]
sq.pl.spatial_scatter(adata_img, color='leiden', shape=None, size=15)

# === 9. SAVE RESULTS ===
print('Saving results...')
adata.write(output_dir / 'imc_analysis.h5ad')

proportions = adata.obs.groupby(['image_id', 'leiden']).size().unstack(fill_value=0)
proportions = proportions.div(proportions.sum(axis=1), axis=0)
proportions.to_csv(output_dir / 'cluster_proportions.csv')

print('Analysis complete!')

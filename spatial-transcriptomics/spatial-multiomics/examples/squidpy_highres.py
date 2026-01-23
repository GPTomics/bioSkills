'''High-resolution spatial analysis with Squidpy'''
import squidpy as sq
import scanpy as sc
import numpy as np

# Load high-resolution spatial data
# Slide-seq, Stereo-seq, or Visium HD
adata = sc.read_h5ad('spatial_highres.h5ad')

# Check spatial coordinates
print(f'Spots/beads: {adata.n_obs}')
print(f'Coordinate range: {adata.obsm["spatial"].min(axis=0)} - {adata.obsm["spatial"].max(axis=0)}')

# Standard preprocessing
sc.pp.filter_cells(adata, min_genes=50)
sc.pp.filter_genes(adata, min_cells=10)
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

# For very high density, bin spots
# bin_size: in coordinate units (usually microns)
# Reduces noise and improves computational efficiency
def bin_spatial(adata, bin_size=50):
    coords = adata.obsm['spatial']
    bins_x = np.floor(coords[:, 0] / bin_size).astype(int)
    bins_y = np.floor(coords[:, 1] / bin_size).astype(int)
    adata.obs['bin_id'] = [f'{x}_{y}' for x, y in zip(bins_x, bins_y)]
    return adata

# adata = bin_spatial(adata, bin_size=50)

# Spatial neighbors
# n_neighs: more neighbors for dense data
sq.gr.spatial_neighbors(
    adata,
    coord_type='generic',
    n_neighs=15,  # Adjust based on spot density
    spatial_key='spatial'
)

# Spatial autocorrelation (Moran's I)
# Identifies spatially variable genes
sq.gr.spatial_autocorr(
    adata,
    mode='moran',
    genes=adata.var_names[:500],  # Top genes
    n_perms=100
)

# Sort by spatial autocorrelation
moran_results = adata.uns['moranI'].sort_values('I', ascending=False)
top_spatial_genes = moran_results.head(50).index.tolist()
print(f'Top spatially variable genes: {top_spatial_genes[:10]}')

# Clustering
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.leiden(adata, resolution=0.5)

# Spatial visualization
sq.pl.spatial_scatter(
    adata,
    color='leiden',
    size=5,  # Small for high density
    save='spatial_clusters.png'
)

# Neighborhood enrichment
sq.gr.nhood_enrichment(adata, cluster_key='leiden')
sq.pl.nhood_enrichment(adata, cluster_key='leiden', save='nhood_enrichment.png')

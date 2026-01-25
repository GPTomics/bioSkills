---
name: bio-imaging-mass-cytometry-spatial-analysis
description: Spatial analysis of cell neighborhoods and interactions in IMC data. Covers neighbor graphs, spatial statistics, and interaction testing. Use when analyzing spatial relationships between cell types, testing for neighborhood enrichment, or identifying cell-cell interaction patterns in imaging mass cytometry data.
tool_type: python
primary_tool: squidpy
---

# Spatial Analysis for IMC

## Build Spatial Graph

```python
import squidpy as sq
import anndata as ad

# Load phenotyped data
adata = ad.read_h5ad('imc_phenotyped.h5ad')

# Ensure spatial coordinates are set
# adata.obsm['spatial'] should contain (x, y) coordinates

# Build spatial neighbor graph
sq.gr.spatial_neighbors(adata, coord_type='generic', delaunay=True)
# Or by distance
sq.gr.spatial_neighbors(adata, coord_type='generic', radius=50)  # 50 pixels

print(f'Built graph with {adata.obsp["spatial_connectivities"].nnz} edges')
```

## Neighborhood Enrichment

```python
# Test if cell types are enriched near each other
sq.gr.nhood_enrichment(adata, cluster_key='cell_type')

# Visualize
sq.pl.nhood_enrichment(adata, cluster_key='cell_type', save='nhood_enrichment.png')

# Get z-scores
zscore = adata.uns['cell_type_nhood_enrichment']['zscore']
# Positive: enriched, Negative: depleted
```

## Co-occurrence Analysis

```python
# Analyze co-occurrence of cell types at multiple distances
sq.gr.co_occurrence(adata, cluster_key='cell_type')

# Plot
sq.pl.co_occurrence(adata, cluster_key='cell_type', save='co_occurrence.png')
```

## Ripley's Statistics

```python
# Ripley's L function for spatial clustering
sq.gr.ripley(adata, cluster_key='cell_type', mode='L')

# Plot
sq.pl.ripley(adata, cluster_key='cell_type', save='ripley.png')

# Interpretation:
# L(r) > r: clustering at distance r
# L(r) < r: dispersion at distance r
# L(r) = r: random distribution
```

## Cell-Cell Interaction

```python
# Permutation test for interactions
sq.gr.interaction_matrix(adata, cluster_key='cell_type', normalized=True)

# Get interaction matrix
interaction = adata.uns['cell_type_interactions']
```

## Custom Neighborhood Analysis

```python
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix

def neighborhood_composition(adata, cluster_key='cell_type'):
    '''Calculate cell type composition of each cell's neighborhood'''

    # Get connectivity matrix
    conn = adata.obsp['spatial_connectivities']
    cell_types = adata.obs[cluster_key]
    type_categories = cell_types.cat.categories

    # One-hot encode cell types
    type_onehot = pd.get_dummies(cell_types).values

    # Neighborhood composition = connectivity * one-hot
    nhood_composition = conn @ type_onehot

    # Normalize to fractions
    nhood_sum = np.array(nhood_composition.sum(axis=1)).flatten()
    nhood_sum[nhood_sum == 0] = 1  # Avoid division by zero
    nhood_frac = nhood_composition / nhood_sum[:, np.newaxis]

    # Add to adata
    for i, ct in enumerate(type_categories):
        adata.obs[f'nhood_frac_{ct}'] = nhood_frac[:, i]

    return nhood_frac

nhood_frac = neighborhood_composition(adata)
```

## Spatial Clustering

```python
# Leiden clustering on spatial + expression
# Weight spatial vs molecular information

# Combined graph
sq.gr.spatial_neighbors(adata, coord_type='generic', radius=30)

# Run spatial Leiden
sc.tl.leiden(adata, adjacency=adata.obsp['spatial_connectivities'],
             resolution=0.5, key_added='spatial_cluster')
```

## Interaction Hotspots

```python
def find_interaction_hotspots(adata, type1, type2, cluster_key='cell_type', radius=50):
    '''Find regions with high interaction between two cell types'''

    # Get cells of each type
    mask1 = adata.obs[cluster_key] == type1
    mask2 = adata.obs[cluster_key] == type2

    spatial = adata.obsm['spatial']

    # For each type1 cell, count nearby type2 cells
    from scipy.spatial import cKDTree

    tree2 = cKDTree(spatial[mask2])

    interaction_scores = np.zeros(mask1.sum())
    for i, (x, y) in enumerate(spatial[mask1]):
        neighbors = tree2.query_ball_point([x, y], r=radius)
        interaction_scores[i] = len(neighbors)

    return interaction_scores

cd8_tumor_interactions = find_interaction_hotspots(adata, 'CD8 T cell', 'Tumor', radius=30)
```

## Visualize Spatial Patterns

```python
import matplotlib.pyplot as plt

# Spatial plot by cell type
sq.pl.spatial_scatter(adata, color='cell_type', size=3, save='spatial_celltypes.png')

# Multiple markers
sq.pl.spatial_scatter(adata, color=['CD8', 'CD4', 'CD68'], size=2, save='spatial_markers.png')

# Highlight specific interaction
fig, ax = plt.subplots(figsize=(10, 10))
spatial = adata.obsm['spatial']

# Background: all cells gray
ax.scatter(spatial[:, 0], spatial[:, 1], c='lightgray', s=1, alpha=0.5)

# Highlight: CD8 and Tumor
for ct, color in [('CD8 T cell', 'red'), ('Tumor', 'blue')]:
    mask = adata.obs['cell_type'] == ct
    ax.scatter(spatial[mask, 0], spatial[mask, 1], c=color, s=5, label=ct)

ax.legend()
ax.set_aspect('equal')
plt.savefig('cd8_tumor_spatial.png', dpi=150)
```

## Statistical Testing

```python
from scipy import stats

def spatial_association_test(adata, type1, type2, cluster_key='cell_type', n_perm=1000):
    '''Permutation test for spatial association between cell types'''

    # Observed interaction count
    sq.gr.nhood_enrichment(adata, cluster_key=cluster_key)
    obs_zscore = adata.uns[f'{cluster_key}_nhood_enrichment']['zscore']

    idx1 = list(adata.obs[cluster_key].cat.categories).index(type1)
    idx2 = list(adata.obs[cluster_key].cat.categories).index(type2)

    observed = obs_zscore[idx1, idx2]

    # The z-score is already normalized, so we can use it directly
    # p-value from z-score
    pvalue = 2 * (1 - stats.norm.cdf(abs(observed)))

    return {'zscore': observed, 'pvalue': pvalue}

result = spatial_association_test(adata, 'CD8 T cell', 'Tumor')
print(f"CD8-Tumor association: z={result['zscore']:.2f}, p={result['pvalue']:.4f}")
```

## Export Results

```python
# Save spatial analysis results
adata.write('imc_spatial_analyzed.h5ad')

# Export neighborhood enrichment
nhood_df = pd.DataFrame(
    adata.uns['cell_type_nhood_enrichment']['zscore'],
    index=adata.obs['cell_type'].cat.categories,
    columns=adata.obs['cell_type'].cat.categories
)
nhood_df.to_csv('neighborhood_enrichment.csv')
```

## Related Skills

- phenotyping - Assign cell types first
- spatial-transcriptomics/spatial-statistics - Similar spatial methods
- single-cell/cell-communication - Interaction concepts

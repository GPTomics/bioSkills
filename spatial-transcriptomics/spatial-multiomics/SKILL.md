---
name: bio-spatial-transcriptomics-spatial-multiomics
description: Analyze high-resolution spatial platforms like Slide-seq, Stereo-seq, and Visium HD. Use when working with subcellular resolution or high-density spatial data.
tool_type: python
primary_tool: squidpy
---

# Spatial Multi-omics Analysis

## Platform Comparison

| Platform | Resolution | Spots/Beads | Coverage |
|----------|------------|-------------|----------|
| Visium | 55 µm | ~5,000 | Tissue-wide |
| Visium HD | 2 µm | ~11M | Subcellular |
| Slide-seq | 10 µm | ~100,000 | High-density |
| Stereo-seq | 0.5 µm | >200M | Subcellular |
| MERFISH | Single-molecule | N/A | Targeted genes |

## Squidpy for High-Resolution Data

```python
import squidpy as sq
import scanpy as sc

# Load spatial data
adata = sc.read_h5ad('spatial_multiomics.h5ad')

# Spatial neighbors (for high-resolution, adjust n_neighs based on density)
sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=10, spatial_key='spatial')

# Spatial autocorrelation (Moran's I)
sq.gr.spatial_autocorr(adata, mode='moran', genes=adata.var_names[:100])

# Neighborhood enrichment analysis
sq.gr.nhood_enrichment(adata, cluster_key='cell_type')
sq.pl.nhood_enrichment(adata, cluster_key='cell_type')

# Ligand-receptor analysis
sq.gr.ligrec(adata, n_perms=100, cluster_key='cell_type')
```

## SpatialData Framework

```python
import spatialdata as sd
from spatialdata_io import read_visium, read_xenium

# Read Visium data
sdata = read_visium('visium_output/')

# Read Xenium data (10x Genomics subcellular)
sdata = read_xenium('xenium_output/')

# Read from Zarr
sdata = sd.read_zarr('experiment.zarr')

# Access different elements
images = sdata.images['morphology']
points = sdata.points['transcripts']
shapes = sdata.shapes['cell_boundaries']
table = sdata.tables['adata']

# Query by region
from spatialdata import bounding_box_query
roi = bounding_box_query(sdata, min_coordinate=[0, 0], max_coordinate=[1000, 1000], axes=['x', 'y'])
```

## Slide-seq/Stereo-seq Processing

```python
# For high-density data, bin spots into hexagonal grids
import numpy as np

# Create hexagonal bins
def hexbin_data(adata, gridsize=50):
    coords = adata.obsm['spatial']
    from matplotlib.pyplot import hexbin
    hb = hexbin(coords[:, 0], coords[:, 1], C=None, gridsize=gridsize, reduce_C_function=np.sum)
    return hb

# Squidpy visualization with hex binning
sq.pl.spatial_scatter(adata, shape='hex', size=50, color='cluster')

# Grid-based spatial neighbors for regular patterns
sq.gr.spatial_neighbors(adata, coord_type='grid', n_rings=1)
```

## Subcellular Analysis (MERFISH/Xenium)

```python
# Transcript-level analysis
# Assign transcripts to compartments
sq.gr.co_occurrence(adata, cluster_key='compartment', spatial_key='spatial')

# Cell segmentation integration
from cellpose import models
model = models.Cellpose(model_type='cyto2')
masks, flows, styles, diams = model.eval(image, diameter=30, channels=[0, 0])

# Map transcripts to cells
def assign_transcripts_to_cells(transcripts_df, masks):
    x, y = transcripts_df['x'].values.astype(int), transcripts_df['y'].values.astype(int)
    transcripts_df['cell_id'] = masks[y, x]
    return transcripts_df[transcripts_df['cell_id'] > 0]
```

## Multi-Modal Integration

```python
# Combine spatial transcriptomics with histology
sq.im.process(adata, layer='image', method='smooth', sigma=2)
sq.im.segment(adata, layer='image', method='watershed', thresh=0.1)

# Extract image features
sq.im.calculate_image_features(
    adata, layer='image', features=['texture', 'summary'],
    key_added='img_features', n_jobs=4
)

# Correlate image features with gene expression
from scipy.stats import pearsonr
for gene in ['marker1', 'marker2']:
    r, p = pearsonr(adata.obs['img_feature'], adata[:, gene].X.flatten())
    print(f'{gene}: r={r:.3f}, p={p:.3e}')
```

## Visium HD Specific

```python
# Visium HD produces bin files at multiple resolutions
# Load 8µm binned data (recommended starting point)
adata = sc.read_h5ad('visium_hd_8um.h5ad')

# Downsample to 16µm if needed for initial analysis
# Original 2µm data available for detailed analysis
```

## Quality Metrics

| Metric | Visium | High-Resolution |
|--------|--------|-----------------|
| Genes/spot | >2000 | >500 |
| UMI/spot | >5000 | >1000 |
| Spatial coverage | >80% | >50% |

## Related Skills

- spatial-transcriptomics/spatial-preprocessing - Standard spatial analysis
- single-cell/preprocessing - scRNA-seq concepts
- spatial-transcriptomics/image-analysis - Morphology processing
- single-cell/cell-annotation - Cell type assignment

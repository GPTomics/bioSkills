---
name: bio-spatial-transcriptomics-spatial-multiomics
description: Analyze high-resolution spatial platforms like Slide-seq, Stereo-seq, and Visium HD. Use when working with subcellular resolution or high-density spatial data.
tool_type: python
primary_tool: Squidpy
---

# Spatial Multi-omics Analysis

## Squidpy for High-Resolution Data

```python
import squidpy as sq
import scanpy as sc

# Load spatial data
adata = sc.read_h5ad('spatial_multiomics.h5ad')

# Spatial neighbors (for high-resolution)
sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=10)

# Spatial autocorrelation
sq.gr.spatial_autocorr(adata, mode='moran', genes=adata.var_names[:100])
```

## SpatialData Framework

```python
import spatialdata as sd

# Read multiple modalities
sdata = sd.read_zarr('experiment.zarr')

# Access different elements
images = sdata.images['morphology']
points = sdata.points['transcripts']
shapes = sdata.shapes['cell_boundaries']
```

## Slide-seq/Stereo-seq Processing

```python
# Bin-based analysis for high-density data
# Aggregate beads/spots into hexagonal bins
sq.pl.spatial_scatter(
    adata,
    shape='hex',
    size=50,  # Bin size in microns
    color='cluster'
)

# Cell type deconvolution per bin
sq.gr.spatial_neighbors(adata, coord_type='grid')
```

## Subcellular Analysis

```python
# Transcript-level analysis
# Count transcripts per cell compartment
sq.gr.co_occurrence(
    adata,
    cluster_key='compartment',  # nucleus, cytoplasm, membrane
    spatial_key='spatial'
)
```

## Related Skills

- **spatial-transcriptomics/basic-analysis** - Standard spatial analysis
- **single-cell/basic-analysis** - scRNA-seq concepts
- **image-analysis** - Morphology processing

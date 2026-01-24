---
name: bio-spatial-transcriptomics-spatial-data-io
description: Load spatial transcriptomics data from Visium, Xenium, MERFISH, Slide-seq, and other platforms using Squidpy and SpatialData. Read Space Ranger outputs, convert formats, and access spatial coordinates. Use when loading Visium, Xenium, MERFISH, or other spatial data.
tool_type: python
primary_tool: squidpy
---

# Spatial Data I/O

Load and work with spatial transcriptomics data from various platforms.

## Required Imports

```python
import squidpy as sq
import scanpy as sc
import anndata as ad
import spatialdata as sd
import spatialdata_io as sdio
```

## Load 10X Visium Data

```python
# Load Space Ranger output (standard method)
adata = sq.read.visium('path/to/spaceranger/output/')
print(f'Loaded {adata.n_obs} spots, {adata.n_vars} genes')

# Spatial coordinates are in adata.obsm['spatial']
print(f"Spatial coords shape: {adata.obsm['spatial'].shape}")

# Image is in adata.uns['spatial']
library_id = list(adata.uns['spatial'].keys())[0]
print(f'Library ID: {library_id}')
```

## Load Visium with Scanpy

```python
# Alternative using Scanpy directly
adata = sc.read_visium('path/to/spaceranger/output/')

# Access tissue image
img = adata.uns['spatial'][library_id]['images']['hires']
scale_factor = adata.uns['spatial'][library_id]['scalefactors']['tissue_hires_scalef']
```

## Load 10X Xenium Data

```python
# Load Xenium output
adata = sq.read.xenium('path/to/xenium/output/')
print(f'Loaded {adata.n_obs} cells')

# Xenium has single-cell resolution
print(f"Cell coordinates: {adata.obsm['spatial'].shape}")
```

## Load with SpatialData (Recommended for New Projects)

```python
import spatialdata_io as sdio

# Load Visium as SpatialData object
sdata = sdio.visium('path/to/spaceranger/output/')
print(sdata)

# Load Xenium
sdata = sdio.xenium('path/to/xenium/output/')

# Access components
table = sdata.tables['table']  # AnnData with expression
shapes = sdata.shapes  # Spatial shapes (spots, cells)
images = sdata.images  # Tissue images
```

## Load MERFISH Data

```python
# MERFISH (Vizgen MERSCOPE)
sdata = sdio.merscope('path/to/merscope/output/')

# Or as AnnData
adata = sq.read.vizgen('path/to/vizgen/output/', counts_file='cell_by_gene.csv', meta_file='cell_metadata.csv')
```

## Load Slide-seq Data

```python
# Slide-seq / Slide-seqV2
adata = sq.read.slideseq('beads.csv', coordinates_file='coords.csv')
```

## Load Nanostring CosMx

```python
# CosMx spatial molecular imaging
sdata = sdio.cosmx('path/to/cosmx/output/')
```

## Load Stereo-seq Data

```python
# Stereo-seq (BGI)
sdata = sdio.stereoseq('path/to/stereoseq/output/')
```

## Load from H5AD with Spatial Coordinates

```python
# If you have h5ad with spatial already stored
adata = sc.read_h5ad('spatial_data.h5ad')

# Verify spatial data exists
if 'spatial' in adata.obsm:
    print('Has spatial coordinates')
if 'spatial' in adata.uns:
    print('Has image data')
```

## Create Spatial AnnData from Scratch

```python
import numpy as np
import pandas as pd

# Expression matrix
X = np.random.poisson(5, size=(1000, 500))

# Spatial coordinates
spatial_coords = np.random.rand(1000, 2) * 1000  # x, y in pixels

# Create AnnData
adata = ad.AnnData(X)
adata.obs_names = [f'spot_{i}' for i in range(1000)]
adata.var_names = [f'gene_{i}' for i in range(500)]
adata.obsm['spatial'] = spatial_coords

# Add minimal spatial metadata for Squidpy
adata.uns['spatial'] = {
    'library_id': {
        'scalefactors': {'tissue_hires_scalef': 1.0, 'spot_diameter_fullres': 50},
    }
}
```

## Access Spatial Coordinates

```python
# Get coordinates as numpy array
coords = adata.obsm['spatial']
x_coords = coords[:, 0]
y_coords = coords[:, 1]

# Get coordinates as DataFrame
coord_df = pd.DataFrame(adata.obsm['spatial'], index=adata.obs_names, columns=['x', 'y'])
```

## Access Tissue Images

```python
# Get high-resolution image
library_id = list(adata.uns['spatial'].keys())[0]
hires_img = adata.uns['spatial'][library_id]['images']['hires']
lowres_img = adata.uns['spatial'][library_id]['images']['lowres']

# Scale factors
scalef = adata.uns['spatial'][library_id]['scalefactors']
print(f"Hires scale: {scalef['tissue_hires_scalef']}")
print(f"Spot diameter: {scalef['spot_diameter_fullres']}")
```

## Convert Between Formats

```python
# SpatialData to AnnData
sdata = sdio.visium('path/to/data/')
adata = sdata.tables['table'].copy()
adata.obsm['spatial'] = np.array(sdata.shapes['spots'][['x', 'y']])

# Save as h5ad
adata.write_h5ad('spatial_converted.h5ad')

# Save SpatialData
sdata.write('spatial_data.zarr')
```

## Load Multiple Samples

```python
# Load and concatenate multiple Visium samples
samples = ['sample1', 'sample2', 'sample3']
adatas = []

for sample in samples:
    adata = sq.read.visium(f'data/{sample}/')
    adata.obs['sample'] = sample
    adatas.append(adata)

# Concatenate
adata_combined = ad.concat(adatas, label='sample', keys=samples)
print(f'Combined: {adata_combined.n_obs} spots')
```

## Subset by Spatial Region

```python
# Select spots in a rectangular region
x_min, x_max = 1000, 2000
y_min, y_max = 1500, 2500

coords = adata.obsm['spatial']
in_region = (coords[:, 0] >= x_min) & (coords[:, 0] <= x_max) & (coords[:, 1] >= y_min) & (coords[:, 1] <= y_max)

adata_region = adata[in_region].copy()
print(f'Selected {adata_region.n_obs} spots')
```

## Related Skills

- spatial-preprocessing - QC and normalization after loading
- spatial-visualization - Plot spatial data
- single-cell/data-io - Non-spatial scRNA-seq data loading

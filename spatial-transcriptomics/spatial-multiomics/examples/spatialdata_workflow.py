'''Multi-modal spatial analysis with SpatialData'''
import spatialdata as sd
from spatialdata_io import xenium, visium_hd
import scanpy as sc

# Load Xenium data (10x subcellular platform)
sdata = xenium.xenium('/path/to/xenium_output')

# Or Visium HD
# sdata = visium_hd.visium_hd('/path/to/visium_hd_output')

# Explore SpatialData structure
print(sdata)
# Contains:
# - images: H&E, DAPI, etc.
# - labels: cell segmentation masks
# - points: transcript locations
# - shapes: cell boundaries
# - table: expression matrix (AnnData)

# Access expression data
adata = sdata.table

# Access transcript coordinates
if 'transcripts' in sdata.points:
    transcripts = sdata.points['transcripts']
    print(f'Total transcripts: {len(transcripts)}')

# Spatial analysis on expression
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.leiden(adata)

# Update table in SpatialData
sdata.table = adata

# Visualize with napari-spatialdata
# from spatialdata_plot import plot
# plot(sdata)

# Aggregate transcripts by cell
# For subcellular analysis
if 'cell_id' in sdata.points['transcripts'].columns:
    cell_counts = sdata.points['transcripts'].groupby('cell_id').size()
    print(f'Transcripts per cell: {cell_counts.mean():.1f}')

# Save processed data
sdata.write('processed.zarr')

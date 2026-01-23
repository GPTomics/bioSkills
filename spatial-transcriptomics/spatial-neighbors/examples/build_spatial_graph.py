'''Build spatial neighbor graph'''

import squidpy as sq
import scanpy as sc

adata = sc.read_h5ad('preprocessed.h5ad')
print(f'Loaded: {adata.n_obs} spots')

# n_neighs=6: Default for hexagonal Visium spots; use 4 for square grids, adjust for spot density
sq.gr.spatial_neighbors(adata, n_neighs=6, coord_type='generic')

conn = adata.obsp['spatial_connectivities']
dist = adata.obsp['spatial_distances']

print(f'\nSpatial neighbor graph:')
print(f'  Edges: {conn.nnz}')
print(f'  Mean neighbors/spot: {conn.nnz / adata.n_obs:.1f}')
print(f'  Mean distance: {dist.data[dist.data > 0].mean():.1f}')

adata.write_h5ad('with_spatial_graph.h5ad')
print('\nSaved to with_spatial_graph.h5ad')

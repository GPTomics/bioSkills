'''Build spatial neighbor graphs under several adjacency definitions and compare them.

The graph IS the model: every downstream spatial statistic is a function of the
weights matrix W, so this script builds W several ways and shows how the choice
changes mean degree and connectivity. Runs on a Squidpy built-in dataset; writes
no files.
'''
# Reference: squidpy 1.4+, scanpy 1.10+, anndata 0.10+, numpy 1.26+ | Verify API if version differs

import squidpy as sq
import numpy as np

adata = sq.datasets.imc()
print(f'loaded {adata.n_obs} cells from squidpy imc dataset')

# Confirm the coordinate scale: median nearest-neighbor spacing tells microns vs pixels apart
sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=1, key_added='nn1')
spacing = np.median(adata.obsp['nn1_distances'].data)
print(f'median nearest-neighbor spacing: {spacing:.1f} coordinate units')


def degree_summary(W):
    deg = np.asarray((W > 0).sum(axis=1)).ravel()
    return deg.mean(), deg.min(), deg.max(), int((deg == 0).sum())


graphs = {}
for k in (6, 15, 30):
    sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=k, key_added=f'knn{k}')
    graphs[f'knn{k}'] = adata.obsp[f'knn{k}_connectivities']

sq.gr.spatial_neighbors(adata, coord_type='generic', delaunay=True, key_added='delaunay')
graphs['delaunay'] = adata.obsp['delaunay_connectivities']

# Prune Delaunay to a physical max edge length to kill spurious cross-gap neighbors;
# the cutoff is a few times the nearest-neighbor spacing measured above
cutoff = 3.0 * spacing
sq.gr.spatial_neighbors(adata, coord_type='generic', delaunay=True, radius=(0.0, cutoff), key_added='delaunay_pruned')
graphs['delaunay_pruned'] = adata.obsp['delaunay_pruned_connectivities']

print('\ngraph                mean_deg  min  max  isolated')
for name, W in graphs.items():
    mean_deg, mn, mx, iso = degree_summary(W)
    print(f'{name:18s}  {mean_deg:7.1f}  {mn:3d}  {mx:3d}  {iso:6d}')

# A downstream statistic inherits all of this. Run Moran's I under each graph and
# compare rankings -- genes that only rank high under one graph are graph-fragile.
genes = adata.var_names[:30].tolist()
for name in ('knn6', 'knn30', 'delaunay_pruned'):
    adata.obsp['spatial_connectivities'] = graphs[name]
    adata.obsp['spatial_distances'] = adata.obsp[f'{name}_distances']
    sq.gr.spatial_autocorr(adata, mode='moran', genes=genes)
    top = adata.uns['moranI'].head(5).index.tolist()
    print(f"top-5 Moran's I genes under {name}: {top}")

'''Load spatial data and reason about platform class, objects, and coordinate frame.

Uses the Squidpy built-in Visium H&E dataset so the script runs without user data
and writes nothing to disk. Real loads swap in squidpy.read.visium('spaceranger_out/')
or a spatialdata_io reader for imaging platforms.
'''
# Reference: squidpy 1.4+, scanpy 1.10+, anndata 0.10+ | Verify API if version differs

import squidpy as sq

adata = sq.datasets.visium_hne_adata()
print(f'Loaded {adata.n_obs} spots, {adata.n_vars} genes')

coords = adata.obsm['spatial']
print(f'Coordinate array shape: {coords.shape}')
print(f'X range: {coords[:, 0].min():.0f} - {coords[:, 0].max():.0f}  (full-res pixels for Visium)')
print(f'Y range: {coords[:, 1].min():.0f} - {coords[:, 1].max():.0f}')

library_id = list(adata.uns['spatial'].keys())[0]
scalef = adata.uns['spatial'][library_id]['scalefactors']
print(f'Library ID: {library_id}')
print(f"Spot diameter: {scalef['spot_diameter_fullres']:.1f} px  (a spot is a 1-10-cell MIXTURE, not a cell)")
print(f"Hires scale factor: {scalef['tissue_hires_scalef']:.4f}  (multiply pixel coords to index the hires image)")

# Visium is a sequencing/capture platform: there is NO per-transcript molecule table.
has_points = 'points' in dir(adata)
print(f'Molecule table present: {has_points}  (capture platforms have none -- nothing to re-segment)')

print(f'Mean total counts/spot: {adata.X.sum(axis=1).mean():.0f}')
print(f'Mean genes detected/spot: {(adata.X > 0).sum(axis=1).mean():.0f}')

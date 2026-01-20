'''Load 10X Visium spatial transcriptomics data'''

import squidpy as sq
import scanpy as sc

adata = sq.read.visium('spaceranger_output/')
print(f'Loaded {adata.n_obs} spots, {adata.n_vars} genes')

coords = adata.obsm['spatial']
print(f'Spatial coordinates shape: {coords.shape}')
print(f'X range: {coords[:, 0].min():.0f} - {coords[:, 0].max():.0f}')
print(f'Y range: {coords[:, 1].min():.0f} - {coords[:, 1].max():.0f}')

library_id = list(adata.uns['spatial'].keys())[0]
print(f'\nLibrary ID: {library_id}')

scalef = adata.uns['spatial'][library_id]['scalefactors']
print(f"Spot diameter: {scalef['spot_diameter_fullres']:.1f} pixels")
print(f"Hires scale factor: {scalef['tissue_hires_scalef']:.4f}")

print(f"\nTotal counts per spot: {adata.X.sum(axis=1).mean():.0f} (mean)")
print(f"Genes detected per spot: {(adata.X > 0).sum(axis=1).mean():.0f} (mean)")

adata.write_h5ad('visium_loaded.h5ad')
print('\nSaved to visium_loaded.h5ad')

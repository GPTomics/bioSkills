'''Extract image features from spatial data'''

import squidpy as sq
import scanpy as sc

adata = sc.read_h5ad('preprocessed.h5ad')
print(f'Loaded: {adata.n_obs} spots')

library_id = list(adata.uns['spatial'].keys())[0]
hires = adata.uns['spatial'][library_id]['images']['hires']
print(f'Image shape: {hires.shape}')

img = sq.im.ImageContainer(hires)

print('\nExtracting image features...')
sq.im.calculate_image_features(
    adata, img,
    features=['summary', 'texture'],
    key_added='img_features',
    spot_scale=1.0,
    n_jobs=4,
)

features = adata.obsm['img_features']
print(f'Extracted {features.shape[1]} features per spot')
print(f'\nFeature statistics:')
print(f'  Mean: {features.mean():.3f}')
print(f'  Std: {features.std():.3f}')

adata.write_h5ad('with_img_features.h5ad')
print('\nSaved to with_img_features.h5ad')

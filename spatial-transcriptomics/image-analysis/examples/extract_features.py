'''Segment a tissue image and extract per-spot image features with Squidpy.

Runs on Squidpy's built-in Visium H&E crop (downloaded once, then cached); writes no files.
Image features describe the pixel patch under each spot -- they do NOT draw cell boundaries.
Watershed here is a fast classical baseline that over-segments H&E; production cell calling
uses Cellpose/StarDist on a nuclear/membrane image or Baysor/proseg on the molecule table.
'''
# Reference: squidpy 1.7+, scikit-image 0.22+, numpy 1.26+, pandas 2.2+ | Verify API if version differs

import warnings
import numpy as np
import pandas as pd
import squidpy as sq
from skimage.measure import regionprops_table

warnings.filterwarnings('ignore')

img = sq.datasets.visium_hne_image_crop()
adata = sq.datasets.visium_hne_adata_crop()
print(f'Loaded {adata.n_obs} spots; image layers: {list(img)}')

# Per-spot image features: summary (intensity stats) + texture (GLCM). Not segmentation.
sq.im.calculate_image_features(adata, img, layer='image', features=['summary', 'texture'],
                               key_added='img_features', n_jobs=1, show_progress_bar=False)
feats = adata.obsm['img_features']
print(f'Extracted {feats.shape[1]} image features for {feats.shape[0]} spots')

# Classical watershed baseline -- adds a segmentation layer to the container.
sq.im.segment(img, layer='image', method='watershed', channel=0, thresh=0.4)
seg = np.asarray(img['segmented_watershed'].values).squeeze().astype(int)
print(f'Watershed labels (baseline, over-segments H&E): {int(seg.max())}')

# Morphology per mask -- the per-cell shape stats QC relies on to spot over/under-segmentation.
props = pd.DataFrame(regionprops_table(seg, properties=['label', 'area', 'eccentricity', 'solidity']))
if len(props):
    print('cell area pct [5,50,95]:', np.percentile(props['area'], [5, 50, 95]).round(1))
print('Done -- treat any derived cell-by-gene matrix as provisional until contamination QC.')

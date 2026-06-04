'''IMC data preprocessing: count-aware hot-pixel removal and cofactor-1 arcsinh.

Demonstrates the count-statistics-correct operations for IMC. The hot-pixel filter
mirrors steinbock --hpf (signed 8-neighbor difference, replace spikes with the neighbor
maximum) -- NOT a median filter, which would zero out the isolated single-positive
pixels that are frequently real signal in a Poisson/zero-inflated regime.
'''
# Reference: numpy 1.26+, scipy 1.12+, tifffile 2024+ | Verify API if version differs
import numpy as np
from scipy import ndimage
import tifffile

def remove_hot_pixels(img, threshold=50):
    # 8-neighbor maximum; a pixel exceeding every neighbor by > threshold (an absolute
    # COUNT difference, not a universal constant -- raise for high-dynamic-range markers)
    # is a spike and is replaced by the neighbor max. Valleys are never filled, so real
    # bright edges and sparse single-positive pixels survive.
    footprint = np.array([[1, 1, 1], [1, 0, 1], [1, 1, 1]])
    neighbor_max = ndimage.maximum_filter(img, footprint=footprint)
    spikes = img - neighbor_max > threshold
    out = img.copy()
    out[spikes] = neighbor_max[spikes]
    return out

def arcsinh_cofactor1(cell_means):
    # cofactor 1 for IMC single-cell means (Hunter 2024 Cytometry A 105:36); the
    # suspension-CyTOF cofactor 5 over-compresses IMC's lower-count cell means
    return np.arcsinh(cell_means / 1.0)

def zscore_cohort(expr, cohort_mean, cohort_std):
    # z-score against COHORT-WIDE statistics, not per-image, so scales stay comparable
    # across samples; per-image scaling destroys cross-sample comparability
    return (expr - cohort_mean) / np.where(cohort_std == 0, 1.0, cohort_std)

img = tifffile.imread('acquisition.tiff').astype(np.float32)  # (channels, y, x) ion counts
print(f'Image shape: {img.shape}')

cleaned = np.stack([remove_hot_pixels(img[c]) for c in range(img.shape[0])])
tifffile.imwrite('cleaned.tiff', cleaned)

# pixel maps are for localization/QC only; quantification happens on per-cell summed
# counts after segmentation, then arcsinh(cofactor 1) and cohort z-score
print('Hot-pixel-filtered stack written; compensate (NNLS) and segment before quantifying.')

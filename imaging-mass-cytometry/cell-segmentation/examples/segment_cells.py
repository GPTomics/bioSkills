'''Whole-cell IMC segmentation with Mesmer, plus the impossible-co-expression QC monitor.

The monitor is the decision-grade output: the fraction of cells co-expressing a
biologically mutually-exclusive marker pair (e.g. CD3 and CD20) is the headline
under-segmentation / lateral-spillover metric -- far more informative than F1/IoU,
which has no tissue ground truth.
'''
# Reference: deepcell 0.12+ (Mesmer), numpy 1.26+, scikit-image 0.22+, tifffile 2024+ | Verify API if version differs
import numpy as np
from skimage import measure
from skimage.segmentation import expand_labels
import tifffile

img = tifffile.imread('compensated.tiff').astype(np.float32)  # (channels, y, x), spillover-compensated
DNA_IDX = 0
MEMBRANE_IDS = [1, 2, 3]   # broadly-expressed membrane markers covering ALL cell types present

def build_membrane(stack, ids):
    # sum pan-membrane markers so no cell type is dark in channel 2; a cell-type-specific
    # sum systematically under-segments the types lacking a marker
    return stack[ids].sum(axis=0)

def segment_mesmer(stack, mpp=1.0):
    from deepcell.applications import Mesmer
    mem = build_membrane(stack, MEMBRANE_IDS)
    two_channel = np.stack([stack[DNA_IDX], mem], axis=-1)[np.newaxis, ...]  # (1, y, x, 2)
    # image_mpp is the TRUE acquisition resolution; Mesmer rescales to its ~0.5 um training
    # resolution, so a wrong mpp degrades every mask
    return Mesmer().predict(two_channel, image_mpp=mpp, compartment='whole-cell')[0, ..., 0]

def nuclear_fallback(nuclear_masks, distance=3):
    # exclusive expansion (stops at the midline between labels) -- a partition, not free
    # dilation; ~3 px at 1 um. Report the radius; accept the macrophage-undercapture bias.
    expanded = expand_labels(nuclear_masks, distance=distance)
    assert expanded.max() == nuclear_masks.max(), 'expansion must not create/destroy cells'
    return expanded

def impossible_coexpression_rate(stack, masks, marker_a, marker_b, pct=75):
    # fraction of cells positive for both members of a mutually-exclusive pair; positivity
    # is per-cell mean above the cohort's pct-th percentile (a sanity gate, not a real call)
    props = measure.regionprops(masks)
    a = np.array([stack[marker_a][p.coords[:, 0], p.coords[:, 1]].mean() for p in props])
    b = np.array([stack[marker_b][p.coords[:, 0], p.coords[:, 1]].mean() for p in props])
    pos_a = a > np.percentile(a, pct)
    pos_b = b > np.percentile(b, pct)
    return float((pos_a & pos_b).mean())

masks = segment_mesmer(img)
print(f'Segmented {masks.max()} cells')

cd3_idx, cd20_idx = 4, 5
rate = impossible_coexpression_rate(img, masks, cd3_idx, cd20_idx)
print(f'CD3+CD20+ rate: {rate:.3f}  (elevated -> suspect under-segmentation / lateral spillover)')

tifffile.imwrite('cell_masks.tiff', masks.astype(np.uint16))

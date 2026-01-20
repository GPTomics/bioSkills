'''Cell segmentation with Cellpose'''
from cellpose import models
import numpy as np
import tifffile
from skimage import measure

# Load preprocessed image
img = tifffile.imread('processed.tiff')
print(f'Image shape: {img.shape}')

# Channel indices (adjust for your panel)
NUCLEAR_CH = 0  # e.g., DNA1
MEMBRANE_CH = 1  # e.g., CD45

# Initialize Cellpose
model = models.Cellpose(model_type='cyto2', gpu=False)

# Prepare input
nuclear = img[NUCLEAR_CH]
membrane = img[MEMBRANE_CH]
img_input = np.stack([membrane, nuclear])

# Segment
print('Running segmentation...')
masks, flows, styles, diams = model.eval(
    img_input,
    channels=[1, 2],
    diameter=40,  # Typical IMC cell diameter in pixels (~10um at 1um/px). Adjust per tissue type.
    flow_threshold=0.4  # Cellpose default; lower (0.1-0.3) for more cells, higher (0.5-0.8) for stringency.
)

print(f'Segmented {masks.max()} cells')

# Extract basic statistics
props = measure.regionprops(masks)
areas = [p.area for p in props]
print(f'Cell area: mean={np.mean(areas):.1f}, median={np.median(areas):.1f}')

# Save
tifffile.imwrite('cell_masks.tiff', masks.astype(np.uint16))
print('Saved masks to cell_masks.tiff')

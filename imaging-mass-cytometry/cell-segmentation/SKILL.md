---
name: bio-imaging-mass-cytometry-cell-segmentation
description: Cell segmentation from multiplexed tissue images. Covers deep learning (Cellpose, Mesmer) and classical approaches for nuclear and whole-cell segmentation. Use when extracting single-cell data from IMC or MIBI images after preprocessing.
tool_type: python
primary_tool: cellpose
---

# Cell Segmentation for IMC

## Cellpose Segmentation

```python
from cellpose import models, io
import numpy as np
import tifffile

# Load image
img = tifffile.imread('processed.tiff')

# Extract nuclear channel (e.g., DNA1)
nuclear_channel = img[0]  # Adjust index based on panel

# Initialize Cellpose model
model = models.Cellpose(model_type='nuclei', gpu=True)

# Run segmentation
masks, flows, styles, diams = model.eval(
    nuclear_channel,
    diameter=30,  # Average nucleus diameter in pixels
    flow_threshold=0.4,
    cellprob_threshold=0.0
)

# masks contains integer labels for each cell
print(f'Cells segmented: {masks.max()}')
```

## Whole-Cell Segmentation with Cellpose

```python
# Use membrane marker for whole-cell
membrane_channel = img[1]  # e.g., CD45

# Combine nuclear and membrane for cyto model
model = models.Cellpose(model_type='cyto2', gpu=True)

# Create 2-channel input [membrane, nuclear]
img_input = np.stack([membrane_channel, nuclear_channel])

masks, flows, styles, diams = model.eval(
    img_input,
    channels=[1, 2],  # [membrane, nuclear]
    diameter=50,
    flow_threshold=0.4
)
```

## Mesmer (DeepCell)

```python
from deepcell.applications import Mesmer

# Initialize Mesmer
app = Mesmer()

# Prepare input: (batch, H, W, 2) - [nuclear, membrane]
img_input = np.stack([nuclear_channel, membrane_channel], axis=-1)
img_input = np.expand_dims(img_input, axis=0)

# Segment
predictions = app.predict(
    img_input,
    image_mpp=1.0,  # Microns per pixel
    compartment='whole-cell'  # or 'nuclear'
)

masks = predictions[0, :, :, 0]
```

## steinbock Segmentation

```bash
# Using steinbock with Cellpose
steinbock segment cellpose \
    --img processed \
    --model cyto2 \
    --channelwise \
    --nuclear-channel 0 \
    --membrane-channel 1 \
    -o masks

# Using steinbock with DeepCell
steinbock segment deepcell \
    --img processed \
    --nuclear-channel 0 \
    --membrane-channel 1 \
    -o masks
```

## Extract Single-Cell Data

```python
from skimage import measure
import pandas as pd

def extract_single_cell_data(img, masks, channel_names):
    '''Extract mean intensity per cell per channel'''

    # Region properties
    props = measure.regionprops(masks)

    # Cell info
    cell_data = []
    intensities = []

    for prop in props:
        # Basic properties
        cell_info = {
            'cell_id': prop.label,
            'area': prop.area,
            'centroid_x': prop.centroid[1],
            'centroid_y': prop.centroid[0],
            'eccentricity': prop.eccentricity
        }
        cell_data.append(cell_info)

        # Mean intensity per channel
        cell_mask = masks == prop.label
        cell_intensities = [img[c][cell_mask].mean() for c in range(len(channel_names))]
        intensities.append(cell_intensities)

    cell_df = pd.DataFrame(cell_data)
    intensity_df = pd.DataFrame(intensities, columns=channel_names)

    return cell_df, intensity_df

cell_info, intensities = extract_single_cell_data(img, masks, channel_names)
print(f'Extracted data for {len(cell_info)} cells')
```

## Quality Control

```python
import matplotlib.pyplot as plt

def qc_segmentation(img, masks, nuclear_channel_idx=0):
    '''Visualize segmentation quality'''

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Nuclear channel
    axes[0].imshow(img[nuclear_channel_idx], cmap='gray')
    axes[0].set_title('Nuclear Channel')

    # Segmentation masks
    axes[1].imshow(masks, cmap='tab20')
    axes[1].set_title(f'Segmentation ({masks.max()} cells)')

    # Overlay
    axes[2].imshow(img[nuclear_channel_idx], cmap='gray')
    axes[2].contour(masks, colors='red', linewidths=0.5)
    axes[2].set_title('Overlay')

    for ax in axes:
        ax.axis('off')

    plt.tight_layout()
    plt.savefig('segmentation_qc.png', dpi=150)
    plt.close()

    # Statistics
    props = measure.regionprops(masks)
    areas = [p.area for p in props]

    print(f'Cells: {len(props)}')
    print(f'Area: mean={np.mean(areas):.1f}, median={np.median(areas):.1f}')

qc_segmentation(img, masks)
```

## Expand Nuclei to Cells

```python
from skimage.segmentation import expand_labels

# If only nuclear segmentation available, expand to approximate cells
nuclear_masks = masks  # From nuclear segmentation
expanded_masks = expand_labels(nuclear_masks, distance=10)

print(f'Expanded masks from nuclei')
```

## Save Results

```python
import tifffile

# Save masks as labeled image
tifffile.imwrite('cell_masks.tiff', masks.astype(np.uint16))

# Save single-cell data
cell_info.to_csv('cell_info.csv', index=False)
intensities.to_csv('cell_intensities.csv', index=False)

# Create combined AnnData
import anndata as ad

adata = ad.AnnData(X=intensities.values)
adata.var_names = channel_names
adata.obs = cell_info

# Add spatial coordinates
adata.obsm['spatial'] = cell_info[['centroid_x', 'centroid_y']].values

adata.write('imc_segmented.h5ad')
```

## Related Skills

- **data-preprocessing** - Prepare images before segmentation
- **phenotyping** - Classify segmented cells
- **spatial-analysis** - Analyze cell spatial relationships

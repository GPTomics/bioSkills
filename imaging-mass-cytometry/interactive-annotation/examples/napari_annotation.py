import numpy as np
import pandas as pd
from pathlib import Path

# Note: napari requires a GUI environment
# This example shows the workflow without actual napari display

# === CONFIGURATION ===
output_dir = Path('annotation_results/')
output_dir.mkdir(exist_ok=True)

# === 1. SIMULATE IMC DATA ===
print('Generating simulated IMC data...')

np.random.seed(42)
height, width = 200, 200
n_channels = 5
channel_names = ['CD45', 'CD3', 'CD68', 'panCK', 'DNA']

# Generate image stack
image_stack = np.random.rand(n_channels, height, width).astype(np.float32)

# Create distinct cell regions
n_cells = 50
segmentation_mask = np.zeros((height, width), dtype=np.int32)

for cell_id in range(1, n_cells + 1):
    cx, cy = np.random.randint(20, 180, 2)
    radius = np.random.randint(5, 15)
    y, x = np.ogrid[:height, :width]
    cell_mask = ((x - cx)**2 + (y - cy)**2) <= radius**2
    segmentation_mask[cell_mask] = cell_id

print(f'Created {n_cells} cells')

# === 2. EXTRACT CELL FEATURES ===
print('\nExtracting cell features...')

cell_features = []
for cell_id in range(1, n_cells + 1):
    mask = segmentation_mask == cell_id
    if not np.any(mask):
        continue

    features = {'cell_id': cell_id}
    for i, name in enumerate(channel_names):
        features[name] = float(np.mean(image_stack[i][mask]))
    features['area'] = int(np.sum(mask))

    cell_features.append(features)

cells_df = pd.DataFrame(cell_features)
print(f'Extracted features for {len(cells_df)} cells')

# === 3. SIMULATE MANUAL ANNOTATIONS ===
print('\nSimulating manual annotations...')

# In practice, this would be done interactively in napari
# Here we simulate by thresholding markers

def assign_cell_type(row):
    if row['CD3'] > 0.5 and row['CD45'] > 0.5:
        return 1  # T cell
    elif row['CD68'] > 0.5:
        return 2  # Macrophage
    elif row['panCK'] > 0.5:
        return 3  # Epithelial
    else:
        return 4  # Other

cells_df['cell_type'] = cells_df.apply(assign_cell_type, axis=1)

cell_type_names = {1: 'T_cell', 2: 'Macrophage', 3: 'Epithelial', 4: 'Other'}
cells_df['cell_type_name'] = cells_df['cell_type'].map(cell_type_names)

print('\nCell type distribution:')
print(cells_df['cell_type_name'].value_counts())

# === 4. CREATE ANNOTATION MASK ===
annotation_mask = np.zeros_like(segmentation_mask)
for _, row in cells_df.iterrows():
    cell_mask = segmentation_mask == row['cell_id']
    annotation_mask[cell_mask] = row['cell_type']

# === 5. VALIDATION SUMMARY ===
print('\n=== ANNOTATION SUMMARY ===')
print(f'Total cells annotated: {len(cells_df)}')
print(f'Cell types: {len(cell_type_names)}')

for ct_id, ct_name in cell_type_names.items():
    count = (cells_df['cell_type'] == ct_id).sum()
    print(f'  {ct_name}: {count} cells')

# === 6. EXPORT ===
cells_df.to_csv(output_dir / 'cell_annotations.csv', index=False)

# Save annotation mask
from PIL import Image
Image.fromarray(annotation_mask.astype(np.uint8)).save(output_dir / 'annotation_mask.png')

print(f'\nResults saved to {output_dir}/')

# === NAPARI EXAMPLE CODE (for reference) ===
napari_code = '''
# To run interactively in napari:

import napari
from skimage import io

# Load data
image_stack = io.imread('imc_image.tiff')
segmentation_mask = io.imread('cell_segmentation.tiff')

# Create viewer
viewer = napari.Viewer()

# Add channels
channel_names = ['CD45', 'CD3', 'CD68', 'panCK', 'DNA']
for i, name in enumerate(channel_names):
    viewer.add_image(image_stack[i], name=name, colormap='gray', blending='additive')

# Add segmentation
viewer.add_labels(segmentation_mask, name='Cells')

# Add annotation layer
import numpy as np
annotation = viewer.add_labels(np.zeros_like(segmentation_mask), name='Cell_Types')

# Annotate by:
# 1. Click cell in 'Cells' layer
# 2. Press number key (1-9) to assign cell type
# 3. Use fill tool to paint annotation

napari.run()
'''

with open(output_dir / 'napari_example.py', 'w') as f:
    f.write(napari_code)

print('Napari example code saved to napari_example.py')

# Cell Segmentation - Usage Guide

## Overview
Cell segmentation identifies individual cells in multiplexed images, providing cell masks essential for downstream single-cell analysis.

## Prerequisites
```bash
pip install cellpose deepcell scikit-image napari
# For GPU acceleration (optional)
pip install cellpose[gpu]
```

## Quick Start
Tell your AI agent what you want to do:
- "Segment cells in my IMC images using Cellpose"
- "Run Mesmer segmentation on my tissue images"
- "Create cell masks from nuclear and membrane channels"

## Example Prompts

### Cellpose Segmentation
> "Segment cells using Cellpose with my DNA channel as nuclear marker"

> "Run Cellpose cyto2 model on my IMC image with diameter 30 pixels"

### Mesmer/DeepCell Segmentation
> "Use Mesmer to segment my tissue image with nuclear and membrane channels"

> "Run DeepCell whole-cell segmentation on my multiplexed image"

### Parameter Tuning
> "Adjust Cellpose flow threshold to reduce oversegmentation"

> "Optimize cell diameter parameter for my tissue type"

### Quality Control
> "Overlay my segmentation mask on the original image for QC"

> "Check my segmentation for oversegmentation and undersegmentation"

## What the Agent Will Do
1. Load preprocessed multichannel image
2. Select appropriate nuclear channel (DNA/Ir-191/193 or histone markers)
3. Optionally select membrane channel (CD45, Na/K-ATPase)
4. Run deep learning model (Cellpose or Mesmer)
5. Post-process masks (remove small objects, fill holes)
6. Save segmentation mask as labeled image

## Tips
- Measure average cell diameter from representative cells for Cellpose
- Cellpose models: nuclei (nuclear only), cyto/cyto2 (whole cell)
- Mesmer is trained specifically on tissue images with better tissue context
- flow_threshold (0.4 default) controls cell separation
- cellprob_threshold (0.0 default) controls detection sensitivity
- Always visually QC segmentation overlaid on original image

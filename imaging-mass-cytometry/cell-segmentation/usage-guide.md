# Cell Segmentation Usage Guide

## Overview

Cell segmentation identifies individual cells in multiplexed images. Critical for downstream single-cell analysis.

## Methods

### Cellpose
- Deep learning, trained on diverse cell types
- Models: nuclei, cyto, cyto2
- Good generalization

### Mesmer (DeepCell)
- Trained specifically on tissue images
- Whole-cell and nuclear options
- Considers tissue context

### Classical
- Watershed, thresholding
- Fast but less accurate
- Good for simple cases

## Channel Selection

### Nuclear Channel
- DNA intercalators (Ir-191/193)
- Histone markers (H3)
- Required for all methods

### Membrane Channel (optional)
- Pan-membrane: CD45, Na/K-ATPase
- Improves whole-cell accuracy

## Parameters

### Diameter
- Cellpose: Average cell diameter in pixels
- Measure from representative cells
- Critical for accuracy

### Thresholds
- flow_threshold: Cell separation (0.4 default)
- cellprob_threshold: Cell probability (0.0 default)

## Quality Control

1. Check oversegmentation (cells split)
2. Check undersegmentation (cells merged)
3. Verify cell areas are reasonable
4. Overlay on image for visual check

## References

- Cellpose: doi:10.1038/s41592-020-01018-x
- Mesmer: doi:10.1038/s41587-021-01094-0

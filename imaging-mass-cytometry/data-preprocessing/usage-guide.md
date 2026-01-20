# Data Preprocessing - Usage Guide

## Overview
IMC/MIBI data requires preprocessing before analysis including hot pixel removal, normalization, and format conversion.

## Prerequisites
```bash
pip install readimc tifffile napari scikit-image
# For steinbock pipeline (Docker-based)
docker pull ghcr.io/bodenmillergroup/steinbock:latest
```

## Quick Start
Tell your AI agent what you want to do:
- "Load MCD files and extract multichannel TIFFs"
- "Remove hot pixels from my IMC images"
- "Normalize channel intensities across samples"

## Example Prompts

### Data Loading
> "Read my MCD file and export each acquisition as OME-TIFF"

> "Convert my Hyperion MCD file to multichannel TIFF images"

### Hot Pixel Removal
> "Detect and remove hot pixels from my IMC images using median filtering"

> "Clean detector artifacts from my MIBI data"

### Normalization
> "Apply percentile normalization to my IMC channels"

> "Arcsinh transform my intensity values for visualization"

### Batch Processing
> "Preprocess all MCD files in my experiment folder"

## What the Agent Will Do
1. Load raw data from MCD/TIFF files using readimc or tifffile
2. Detect hot pixels by comparing to local median filter
3. Replace hot pixels with local median values
4. Apply normalization (percentile, z-score, or min-max)
5. Optionally apply arcsinh transformation for visualization
6. Save preprocessed images as OME-TIFF

## Tips
- MCD files contain all acquisitions; extract each ROI separately
- Hot pixels appear as bright single-pixel spots from detector noise
- Percentile normalization (1st-99th) is robust to outliers
- Arcsinh transform (cofactor 5) is standard for visualization
- steinbock provides a complete Docker-based preprocessing pipeline

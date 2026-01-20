# Data Preprocessing Usage Guide

## Overview

IMC/MIBI data requires preprocessing before analysis. Key steps include hot pixel removal, normalization, and format conversion.

## Data Formats

| Format | Source | Description |
|--------|--------|-------------|
| MCD | Hyperion | Contains all acquisitions |
| TIFF | Export | Multi-channel images |
| OME-TIFF | Standard | Metadata-rich TIFF |

## Preprocessing Steps

1. **Extract images** from MCD files
2. **Remove hot pixels** (detector artifacts)
3. **Normalize** intensity values
4. **Transform** (arcsinh for visualization)

## Hot Pixels

Caused by detector noise. Appear as bright spots in single pixels.
- Median filter comparison detects them
- Replace with local median

## Normalization Options

- **Percentile**: Scale to 1st-99th percentile
- **Z-score**: Mean=0, SD=1 per channel
- **Min-max**: Scale to 0-1

## Tools

- **steinbock**: Complete pipeline (Docker-based)
- **readimc**: Python MCD reading
- **napari**: Visualization
- **tifffile**: TIFF I/O

## References

- steinbock: doi:10.1038/s41592-022-01459-z
- readimc: https://github.com/BodenmillerGroup/readimc

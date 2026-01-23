# Spatial Data I/O - Usage Guide

## Overview

This skill covers loading spatial transcriptomics data from various platforms including 10X Visium, Xenium, MERFISH, Slide-seq, and more using Squidpy and SpatialData.

## Prerequisites

```bash
pip install squidpy spatialdata spatialdata-io scanpy anndata
```

## Quick Start

Tell your AI agent what you want to do:

- "Load my Visium data from the Space Ranger output"
- "Read this Xenium experiment"
- "Load my MERFISH spatial data"

## Example Prompts

### Visium
> "Load my 10X Visium data"

> "Read the Space Ranger output in this folder"

### Xenium
> "Load my Xenium data"

> "Read single-cell spatial data from Xenium"

### Other Platforms
> "Load my MERFISH data"

> "Read Slide-seq bead coordinates"

> "Import CosMx data"

## What the Agent Will Do

1. Identify the appropriate reader for your platform
2. Load expression data into AnnData format
3. Extract spatial coordinates
4. Load tissue images if available
5. Return a properly formatted spatial object

## Tips

- **Visium** uses spots (~55um), while Xenium/MERFISH have single-cell resolution
- **SpatialData** is the newer format that handles multi-modal spatial data well
- **Library ID** is needed to access images in Visium data
- **Coordinates** are stored in `adata.obsm['spatial']` as (x, y) pairs

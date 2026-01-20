# IMC Pipeline Usage Guide

## Overview

This workflow processes imaging mass cytometry data from raw acquisitions through cell segmentation, phenotyping, and spatial analysis.

## When to Use This Pipeline

- Imaging mass cytometry (IMC) tissue analysis
- CODEX/MIBI-TOF multiplexed imaging
- High-plex tissue phenotyping
- Tumor microenvironment studies
- Spatial cell interaction analysis

## Required Inputs

1. **Raw images** - MCD files or OME-TIFF stacks
2. **Panel file** - Channel to marker mapping
3. **Sample metadata** - Conditions, patient IDs

## Panel File Format

```csv
channel,name,full_name,keep,segment
1,DNA1,Ir191,1,1
2,DNA2,Ir193,1,1
3,CD45,CD45,1,0
4,CD3,CD3,1,0
5,CD8,CD8,1,0
```

## Pipeline Steps

### 1. Data Preprocessing
- Extract images from MCD files
- Hot pixel filtering
- Optional: spillover correction

### 2. Cell Segmentation
- Cellpose: Deep learning, good for varied morphologies
- Mesmer: Trained on tissue images, uses nuclear+membrane
- Choose based on tissue type

### 3. Single-cell Quantification
- Mean intensities per cell per marker
- Cell properties (area, eccentricity)
- Neighbor relationships

### 4. Clustering & Phenotyping
- Arcsinh transformation (cofactor 5)
- Leiden/FlowSOM clustering
- Annotate based on marker expression

### 5. Spatial Analysis
- Neighborhood enrichment
- Co-occurrence statistics
- Ligand-receptor interactions

## Segmentation Parameters

| Tissue Type | Method | Diameter | Notes |
|-------------|--------|----------|-------|
| FFPE | Cellpose | 15-20 | Smaller nuclei |
| Frozen | Cellpose | 20-30 | Larger cells |
| Dense tumor | Mesmer | Auto | Better separation |
| Sparse | Cellpose | 25-35 | Larger diameter |

## Quality Metrics

### Segmentation Quality
- Check segmentation masks visually
- Cell size distribution (median 100-500 px)
- Avoid over/under-segmentation

### Expression Quality
- Markers show expected patterns
- No obvious batch effects in UMAP
- Controls express expected markers

## Common Issues

### Poor segmentation
- Adjust diameter parameter
- Try different nuclear channels
- Use membrane for Mesmer

### Hot pixel artifacts
- Lower hot pixel threshold
- Check raw images

### No spatial signal
- Verify neighbor detection distance
- Check if cells are too sparse

### Batch effects
- Use batch-aware methods (Harmony, scVI)
- Include multiple ROIs per condition

## Output Files

| File | Description |
|------|-------------|
| masks/ | Cell segmentation masks |
| intensities/ | Single-cell marker expression |
| imc_analysis.h5ad | Complete analysis object |
| cell_type_proportions.csv | Cell type frequencies |
| umap_celltypes.png | Cluster visualization |
| spatial_celltypes.png | Spatial cell maps |
| neighborhood_enrichment.png | Spatial interactions |

## References

- steinbock: doi:10.1186/s12859-022-04716-7
- Cellpose: doi:10.1038/s41592-020-01018-x
- squidpy: doi:10.1038/s41592-021-01358-2
- imcRtools: Bioconductor package

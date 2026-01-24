# Flow Cytometry Pipeline Usage Guide

## Overview

This workflow processes flow cytometry data from raw FCS files through compensation, transformation, clustering or gating, and differential analysis.

## Prerequisites

```r
BiocManager::install(c('CATALYST', 'flowCore', 'diffcyt'))
install.packages(c('FlowSOM', 'uwot'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Run the cytometry pipeline on my FCS files"
- "Cluster my CyTOF data and find differential populations"
- "Analyze my immunophenotyping experiment"

## Example Prompts

### Basic Analysis
> "I have FCS files from a CyTOF experiment, run the full pipeline"

> "Cluster my flow cytometry data with FlowSOM"

### Differential Analysis
> "Compare immune populations between treatment and control groups"

> "Find differentially abundant cell types in my flow data"

### Customization
> "Use manual gating instead of clustering for my data"

> "Run differential state analysis on my T cell markers"

## When to Use This Pipeline

- Multi-color flow cytometry panels
- CyTOF/mass cytometry experiments
- Immunophenotyping studies
- Treatment response analysis
- Biomarker discovery

## Required Inputs

1. **FCS files** - One per sample
2. **Panel file** - Channel to marker mapping with marker classes
3. **Sample metadata** - Conditions, patient IDs, batches

## Panel File Format

```csv
fcs_colname,antigen,marker_class
FSC-A,FSC,none
SSC-A,SSC,none
CD45,CD45,type
CD3,CD3,type
CD4,CD4,type
Ki67,Ki67,state
IFNg,IFNg,state
```

**Marker classes:**
- `type`: Lineage markers for clustering/phenotyping
- `state`: Functional markers for differential state analysis
- `none`: Scatter or non-analyzed channels

## Pipeline Steps

### 1. Data Loading
- Read FCS files into R
- Validate panel and metadata

### 2. Compensation
- Apply spillover matrix
- Essential for fluorescence-based flow
- Not needed for CyTOF

### 3. Transformation
- Arcsinh transformation (standard)
- Cofactor: 5 for CyTOF, 150 for flow
- Normalizes distributions

### 4. QC
- Check cell counts per sample
- Visualize expression distributions
- Identify outlier samples

### 5. Clustering or Gating
- **Clustering**: FlowSOM, Phenograph (unsupervised)
- **Gating**: Manual population identification

### 6. Dimensionality Reduction
- UMAP for visualization
- Based on lineage markers

### 7. Differential Analysis
- **DA**: Are population proportions different?
- **DS**: Is marker expression different within populations?

## Analysis Approaches

### Approach 1: Unsupervised Clustering
Best for: High-dimensional panels, discovery

1. Cluster on lineage markers
2. Annotate clusters based on expression
3. Test differential abundance

### Approach 2: Manual Gating
Best for: Well-defined populations, lower dimensions

1. Gate populations hierarchically
2. Export gated populations
3. Compare frequencies

## Cofactor Guidelines

| Platform | Typical Cofactor |
|----------|------------------|
| Flow cytometry | 150-500 |
| CyTOF | 5 |
| Spectral flow | 500-1000 |

## Statistical Considerations

### Differential Abundance
- Need biological replicates (n ≥ 3)
- Uses edgeR (count-based)
- Reports fold change in frequency

### Differential State
- Tests marker expression within clusters
- Uses limma
- Adjust for multiple testing (clusters × markers)

## Common Issues

### Compensation problems
- Check single-stain controls
- Verify spillover matrix
- Look for negative values

### Poor clustering
- Adjust number of clusters (K)
- Check if markers are appropriate
- Try different resolution

### No significant results
- Check sample size (n ≥ 3)
- Verify biological effect expected
- Adjust thresholds carefully

## Output Files

| File | Description |
|------|-------------|
| cytometry_analysis.rds | Complete analysis object |
| da_results.csv | Differential abundance results |
| ds_results.csv | Differential state results |
| umap_clusters.png | Cluster visualization |
| abundance_boxplots.png | Population frequencies |
| da_volcano.png | Differential abundance plot |

## Tips

- **Panel file**: Essential for correct channel-to-marker mapping
- **Cofactor**: Use 5 for CyTOF, 150 for conventional flow
- **Batch effects**: Include batch in metadata if samples were processed separately
- **Replicates**: Minimum 3 biological replicates per condition for statistical testing
- **Marker classes**: Separate lineage (type) from functional (state) markers

## References

- CATALYST: doi:10.1101/218826
- diffcyt: doi:10.1038/s41467-017-00707-4
- FlowSOM: doi:10.1002/cyto.a.22625

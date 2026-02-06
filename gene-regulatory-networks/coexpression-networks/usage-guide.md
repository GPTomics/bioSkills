# Co-expression Networks - Usage Guide

## Overview

Build weighted gene co-expression networks to identify modules of co-regulated genes and relate them to phenotypes. WGCNA constructs a scale-free network from expression data, groups genes into modules based on co-expression patterns, and correlates module eigengenes with sample traits to find biologically relevant modules and hub genes.

## Prerequisites

```r
# R packages
install.packages('WGCNA')
install.packages('CEMiTool')
BiocManager::install('hdWGCNA')

# Python alternative
pip install PyWGCNA
```

Minimum 20 samples recommended for reliable module detection (absolute floor of 15).

## Quick Start

Tell your AI agent what you want to do:
- "Build a co-expression network from my RNA-seq data"
- "Find gene modules correlated with treatment response"
- "Identify hub genes in the most significant module"
- "Run WGCNA on my bulk RNA-seq count matrix"
- "Find co-expression modules in my single-cell data"

## Example Prompts

### Module Detection
> "I have normalized RNA-seq counts from 30 samples. Build a WGCNA network and find gene modules."

> "Run CEMiTool on my expression data for an automated co-expression analysis."

### Module-Trait Relationships
> "Correlate my WGCNA modules with survival time and treatment group."

> "Which modules are significantly associated with disease status?"

### Hub Genes
> "Find hub genes in the module most correlated with my phenotype."

> "Export the turquoise module network for visualization in Cytoscape."

### Single-Cell
> "Run hdWGCNA on my Seurat object to find co-expression modules per cell type."

## What the Agent Will Do

1. Load and filter expression data for low-variance genes
2. Select soft-thresholding power based on scale-free topology fit (R^2 > 0.85)
3. Construct network and detect gene modules
4. Calculate module eigengenes and correlate with sample traits
5. Identify hub genes based on module membership and trait significance
6. Visualize module-trait heatmap and export networks

## Tips

- **Sample size matters** - WGCNA needs at least 15 samples, 20+ recommended for robust results
- **Soft power selection** - Choose the lowest power where scale-free R^2 exceeds 0.85; if it never reaches 0.85, aim for the plateau above 0.80
- **Batch effects** - Remove batch effects before network construction; batch-driven modules are artifacts
- **Gene filtering** - Use the top 5000-10000 most variable genes to reduce noise and speed computation
- **CEMiTool for quick analysis** - Automates soft-threshold selection and module detection; good for exploratory analysis
- **hdWGCNA for single-cell** - Creates metacells to address scRNA-seq sparsity before applying WGCNA

## Related Skills

- scenic-regulons - TF regulon inference from scRNA-seq with pySCENIC
- differential-networks - Compare networks between conditions with DiffCorr
- differential-expression/deseq2-basics - DE analysis for gene prioritization
- differential-expression/batch-correction - Remove batch effects before network construction
- pathway-analysis/go-enrichment - Functional enrichment of gene modules

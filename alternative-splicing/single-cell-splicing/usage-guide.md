# Single-Cell Splicing - Usage Guide

## Overview
Analyze alternative splicing at single-cell resolution. Handles the sparsity of scRNA-seq data using probabilistic models and identifies cell-type-specific splicing patterns.

## Prerequisites
```bash
# BRIE2
pip install brie

# Additional dependencies
pip install scanpy anndata pandas numpy

# leafcutter2 (for intron clustering approach)
# git clone https://github.com/davidaknowles/leafcutter
```

## Quick Start
Tell your AI agent what you want to do:
- "Analyze splicing heterogeneity in my scRNA-seq data"
- "Find cell-type-specific splicing patterns"
- "Calculate PSI values for single cells using BRIE2"
- "Identify differential splicing between cell populations"

## Example Prompts

### BRIE2 Analysis
> "Run BRIE2 to estimate PSI values for splicing events in my 10X scRNA-seq data."

> "Use BRIE2 to find variable splicing events across my single cells."

### Cell-Type Specific
> "Find splicing events that differ between cell types in my annotated scRNA-seq."

> "Calculate mean PSI per cell type and identify variable events."

### leafcutter2
> "Run leafcutter2 for intron cluster analysis on my single-cell BAM."

> "Detect NMD-inducing splicing events in my scRNA-seq data."

## What the Agent Will Do
1. Prepare splicing event annotations from GTF
2. Count junction reads per cell from BAM files
3. Run BRIE2 probabilistic PSI estimation
4. Aggregate by cell type for comparison
5. Identify variable splicing events across populations

## Tips
- BRIE2 handles scRNA-seq sparsity with probabilistic modeling
- Consider pseudobulk approach for better statistical power
- Avoid Whippet.jl (incompatible with Julia 1.9+)
- leafcutter2 (April 2025) adds NMD detection capability
- 10X 3' data has bias against 5' splicing events; consider 5' kit
- Filter cells for quality before splicing analysis
- Require at least 50 cells with reads per event for reliability

## Related Skills
- single-cell/preprocessing - QC before splicing analysis
- single-cell/clustering - Cell type annotation
- splicing-quantification - Bulk RNA-seq comparison

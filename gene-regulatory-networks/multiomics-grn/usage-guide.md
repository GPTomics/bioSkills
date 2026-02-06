# Multiomics GRN Inference - Usage Guide

## Overview

Build enhancer-driven gene regulatory networks by integrating single-cell RNA-seq and ATAC-seq data. SCENIC+ extends the SCENIC framework to discover eRegulons: triplets linking transcription factors to their target enhancers and downstream regulated genes. This enables identification of cis-regulatory programs driving cell identity.

## Prerequisites

```bash
# SCENIC+ and dependencies
pip install scenicplus pycisTopic

# Peak calling
pip install macs3

# Additional
pip install scanpy loompy matplotlib seaborn

# Download cisTarget databases for SCENIC+
# See: https://resources.aertslab.org/cistarget/
```

```r
# FigR alternative (R)
devtools::install_github('buenrostrolab/FigR')
BiocManager::install(c('Signac', 'Seurat'))
```

## Quick Start

Tell your AI agent what you want to do:
- "Infer gene regulatory networks from my 10x multiome data"
- "Run SCENIC+ on my paired scRNA+scATAC dataset"
- "Find eRegulons linking TFs to enhancers and target genes"
- "Build a multiomics GRN from my CellRanger ARC output"

## Example Prompts

### SCENIC+ Pipeline
> "I have 10x multiome data processed with CellRanger ARC. Run SCENIC+ to identify eRegulons."

> "Infer enhancer-driven regulatory networks from my paired scRNA-seq and scATAC-seq data."

### Interpretation
> "Which TFs have the most enhancer targets in my eRegulon results?"

> "Show activating vs repressive eRegulons for each cell type."

### FigR Alternative
> "Run FigR on my Seurat multiome object to find TF-gene regulatory links."

## What the Agent Will Do

1. Load paired scRNA-seq and scATAC-seq data
2. Call peaks from ATAC fragments with MACS3
3. Run cisTopic for topic modeling on scATAC regions
4. Link accessible regions to target genes (region-to-gene)
5. Connect TFs to enhancers via motif enrichment (TF-to-region)
6. Assemble eRegulons (TF-enhancer-gene triplets)
7. Score eRegulon activity per cell and visualize

## Tips

- **Memory requirements** - SCENIC+ needs significant RAM: 64 GB for 20K cells, 128 GB+ for 50K+ cells
- **cisTopic topics** - Test multiple topic numbers (10-50) and select by coherence metrics
- **Peak calling** - Use MACS3 with BEDPE format on the fragments file for the region universe
- **Activating vs repressive** - eRegulons with (+) suffix have positive TF-target correlation; (-) suffix indicates repressive regulation
- **FigR is lighter weight** - If SCENIC+ is too resource-intensive, FigR provides a simpler alternative for TF-gene regulatory inference from multiome data

## Related Skills

- scenic-regulons - RNA-only regulon inference with pySCENIC
- single-cell/scatac-analysis - scATAC-seq preprocessing with Signac and ArchR
- atac-seq/atac-peak-calling - Peak calling for chromatin accessibility
- chip-seq/motif-analysis - Motif enrichment for TF identification

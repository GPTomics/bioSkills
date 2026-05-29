# Multiomics GRN Inference - Usage Guide

## Overview

Build enhancer-driven gene regulatory networks (eGRNs) by integrating single-cell RNA-seq and ATAC-seq. SCENIC+ discovers eRegulons: triplets linking transcription factors to target enhancers and downstream genes. The core advance is that scATAC nominates the actual distal enhancers active in the cells, so an edge is a TF -> region -> gene triplet with the region as a validatable anchor. But every inference arrow leaks (motif != binding, accessible != this-TF, peak-gene correlation != control), the peak-to-gene linking step is confounded by cell-type composition, and metacell aggregation introduces pseudo-replication -- so an eGRN is a prioritized hypothesis list, not a wiring diagram. Method choice hinges on paired vs unpaired data and on whether the downstream goal is perturbation simulation (build a CellOracle base GRN here, simulate in perturbation-simulation).

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
> "Run FigR on my Seurat multiome object to find TF-gene regulatory links via DORCs."

### Unpaired or RNA-only
> "I have separate scRNA-seq and scATAC-seq experiments, not paired. Integrate them with GLUE before inferring the GRN."

> "I only have scRNA-seq. Build a GRN using a prebuilt CellOracle base GRN."

## What the Agent Will Do

1. Load paired scRNA-seq and scATAC-seq data
2. Call peaks from ATAC fragments with MACS3
3. Run cisTopic for topic modeling on scATAC regions
4. Link accessible regions to target genes (region-to-gene)
5. Connect TFs to enhancers via motif enrichment (TF-to-region)
6. Assemble eRegulons (TF-enhancer-gene triplets)
7. Score eRegulon activity per cell and visualize

## Tips

- Accessibility defines enhancers, not function - an eGRN edge means "TF motif in an accessible peak linked to a gene," a hypothesis to validate with ChIP/CRISPRi, not established regulation.
- Watch the composition confound - genome-wide peak-gene correlation recovers cell-type co-markers; restrict to the distance window, use a matched null (Signac) or within-type correlation, and require a motif.
- Metacell p-values are ranking scores - KNN-smoothed metacells are not independent; do not report their p-values as calibrated significance.
- Use the current SCENIC+ Snakemake workflow - the manual-object API and pre-2024 tutorials are stale; consensus peak calling needs cell-type labels first.
- Paired vs unpaired drives method choice - paired multiome -> SCENIC+/Pando/TRIPOD/FigR directly; unpaired -> integrate with GLUE first (links are capped by pairing accuracy); RNA-only -> CellOracle prebuilt base GRN.
- Don't over-single-out a TF - motif families (GATA, ETS, bZIP) share motifs; report the family unless orthogonal evidence singles out a member.
- Resource planning - SCENIC+ is RAM-heavy (~64 GB for 20K cells, 128 GB+ beyond); FigR/DIRECT-NET are lighter alternatives.

## Related Skills

- scenic-regulons - RNA-only regulon inference with pySCENIC
- perturbation-simulation - in silico TF perturbation built on a CellOracle base GRN
- single-cell/scatac-analysis - scATAC-seq preprocessing with Signac and ArchR
- atac-seq/atac-peak-calling - peak calling for chromatin accessibility
- chip-seq/motif-analysis - motif enrichment and TF binding for validation

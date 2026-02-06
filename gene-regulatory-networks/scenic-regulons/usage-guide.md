# SCENIC Regulons - Usage Guide

## Overview

Infer transcription factor regulons from single-cell RNA-seq data using the pySCENIC three-step pipeline. Discovers co-expression modules between TFs and candidate targets with GRNBoost2, prunes these to retain only links supported by cis-regulatory motif enrichment (cisTarget), and scores regulon activity per cell with AUCell. Identifies master regulators of cell identity and enables TF-activity-based cell clustering.

## Prerequisites

```bash
# Create dedicated environment (pySCENIC tested on Python 3.10)
conda create -n scenic python=3.10
conda activate scenic
pip install pyscenic loompy scanpy matplotlib seaborn

# Download required databases (human hg38 example)
# Ranking databases (~1.5 GB each)
wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather

# Motif annotations
wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl

# TF list
wget https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt
```

## Quick Start

Tell your AI agent what you want to do:
- "Identify transcription factor regulons in my single-cell data"
- "Run pySCENIC on my scRNA-seq to find master regulators"
- "Score TF activity per cell and find cell-type-specific regulons"
- "Which TFs drive the identity of each cluster?"

## Example Prompts

### Full Pipeline
> "I have a preprocessed scRNA-seq h5ad file. Run the full pySCENIC pipeline to identify TF regulons."

> "Run GRNBoost2, cisTarget pruning, and AUCell scoring on my loom file."

### Regulon Interpretation
> "Find which regulons are specific to each cell type using RSS scores."

> "Binarize regulon activity and show the fraction of cells with active regulons per cluster."

### Visualization
> "Create a heatmap of regulon activity across cell types."

> "Plot the top 3 regulons per cell type on the UMAP embedding."

## What the Agent Will Do

1. Convert expression data to loom format if needed
2. Run GRN inference with GRNBoost2 (using arboreto_with_multiprocessing.py to avoid dask issues)
3. Prune co-expression modules by cis-regulatory motif enrichment with cisTarget
4. Score regulon activity per cell with AUCell
5. Calculate regulon specificity scores per cell type
6. Visualize regulon activity on UMAP and as heatmaps

## Tips

- **Dask compatibility** - Native Arboreto/GRNBoost2 is broken with dask >= 2.0; always use the arboreto_with_multiprocessing.py script bundled with pySCENIC
- **Python version** - Use Python 3.10 in a dedicated conda environment to avoid dependency conflicts
- **Subsample for speed** - Step 1 (GRNBoost2) is the bottleneck; subsample to 5000-10000 cells, then score regulons on the full dataset in Step 3
- **Database matching** - Ranking databases must match genome build (hg38/mm10) and gene naming convention
- **Loom format** - pySCENIC expects loom input; convert from h5ad using loompy
- **Gene filtering** - Remove genes expressed in fewer than 3 cells before running

## Related Skills

- multiomics-grn - Enhancer-driven GRNs from paired scRNA+scATAC with SCENIC+
- coexpression-networks - Bulk co-expression analysis with WGCNA
- single-cell/clustering - Cluster cells before regulon analysis
- single-cell/preprocessing - QC and normalization of scRNA-seq inputs

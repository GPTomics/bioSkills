# SCENIC Regulons - Usage Guide

## Overview

Infer transcription factor regulons from single-cell RNA-seq using the pySCENIC three-step pipeline. GRNBoost2 finds TF-target co-expression, cisTarget prunes it to direct targets by requiring the TF's motif to be enriched in their regulatory regions, and AUCell scores per-cell regulon activity. The key idea: the motif-pruning step is what turns undirected co-expression into directed, motif-supported regulation, and AUCell activity is a rank statistic, so regulon activity is distinct from TF expression (it can be high while the TF transcript is dropout-zero). Regulons are strong directed hypotheses, not proof of causal regulation. For paired scRNA+scATAC use multiomics-grn; for bulk data and TF protein-activity use grn-inference.

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

### Reproducibility and Comparison
> "Run the GRN step with several seeds and keep only the regulon links that recur in most runs."

> "Compare regulon activity between treated and control cells on one integrated SCENIC run."

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

- Motif pruning is the point - step 1 alone is co-expression; only the cisTarget step yields directed, motif-supported regulons. Do not call unpruned modules regulons.
- Activity is not expression - AUCell reports whether the TF's target program is coordinately high-ranked in a cell, robust to dropout; report regulon AUC, not TF mRNA.
- Dask compatibility - native arboreto/GRNBoost2 breaks on dask >= 2.x; use the arboreto_with_multiprocessing.py script bundled with pySCENIC.
- Reproducibility - GRNBoost2 is stochastic; fix the seed, and for confident regulons run the GRN step many times and keep recurring links.
- Database matching - ranking DB, motif2TF annotations, and gene IDs must all match species/assembly/symbol namespace; mismatches give empty regulons. Run ctx with both search-space DBs.
- Compare within one run - raw AUC is population-relative; for cross-condition comparison run SCENIC once on the integrated object and check condition regulons are not batch artifacts.
- QC first - remove doublets and tiny clusters before SCENIC; both inflate spurious regulons.

## Related Skills

- multiomics-grn - enhancer-driven eRegulons from paired scRNA+scATAC with SCENIC+
- grn-inference - bulk GRN inference and VIPER TF protein-activity
- coexpression-networks - undirected co-expression modules (what step 1 produces alone)
- single-cell/clustering - cluster cells before regulon analysis
- single-cell/preprocessing - QC and normalization of scRNA-seq inputs
- single-cell/doublet-detection - remove doublets that inflate spurious regulons

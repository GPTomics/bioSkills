---
name: bio-gene-regulatory-networks-multiomics-grn
description: Build enhancer-driven gene regulatory networks by integrating single-cell RNA-seq and ATAC-seq data using SCENIC+ to identify eRegulons linking transcription factors to enhancers and target genes. Use when analyzing 10x multiome or paired scRNA+scATAC data to infer cis-regulatory GRNs.
tool_type: python
primary_tool: SCENIC+
---

## Version Compatibility

Reference examples tested with: Cell Ranger 8.0+, MACS3 3.0+, matplotlib 3.8+, pandas 2.2+, scanpy 1.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Multiomics GRN Inference

**"Build an enhancer-driven gene regulatory network from my multiome data"** → Integrate scRNA-seq and scATAC-seq to identify eRegulons: transcription factor-enhancer-target gene triplets linking TF binding to chromatin accessibility and gene expression changes.
- Python: SCENIC+ pipeline with `scenicplus` for eRegulon assembly
- Python: `pycisTopic` for topic modeling of scATAC-seq regions

Build enhancer-driven gene regulatory networks from paired single-cell RNA-seq and ATAC-seq data. SCENIC+ extends SCENIC by linking TFs to their enhancers and target genes through eRegulons.

## SCENIC+ Overview

| Component | Tool | Purpose |
|-----------|------|---------|
| Topic modeling | cisTopic | Identify cis-regulatory topics from scATAC |
| Region-to-gene | SCENIC+ | Link accessible regions to target genes |
| TF-to-region | SCENIC+ | Connect TFs to their enhancer targets |
| eRegulon assembly | SCENIC+ | Combine into TF-enhancer-gene triplets |

## Input Preparation

### From 10x Multiome (CellRanger ARC)

```python
import scanpy as sc
import pycisTopic
from pycisTopic.cistopic_class import create_cistopic_object_from_fragments

# Load RNA from CellRanger ARC output
adata_rna = sc.read_10x_h5('filtered_feature_bc_matrix.h5', gex_only=True)
adata_rna.var_names_make_unique()
sc.pp.filter_cells(adata_rna, min_genes=200)
sc.pp.filter_genes(adata_rna, min_cells=3)

# Load ATAC fragments
fragments_file = 'atac_fragments.tsv.gz'
```

### Call Peaks with MACS3

```python
import subprocess

# Call peaks from fragments file for region universe
subprocess.run([
    'macs3', 'callpeak',
    '-t', 'atac_fragments.tsv.gz',
    '-f', 'BEDPE',
    '--nomodel', '--shift', '-75', '--extsize', '150',
    '-g', 'hs',
    '--keep-dup', 'all',
    '-n', 'multiome_peaks',
    '--outdir', 'macs3_output'
], check=True)
```

### Create cisTopic Object

```python
from pycisTopic.cistopic_class import create_cistopic_object_from_fragments
from pycisTopic.lda_models import run_cgs_models

# Create binary accessibility matrix from fragments
path_to_regions = 'macs3_output/multiome_peaks_peaks.narrowPeak'
cistopic_obj = create_cistopic_object_from_fragments(
    path_to_fragments=fragments_file,
    path_to_regions=path_to_regions,
    path_to_blacklist='hg38-blacklist.v2.bed',
    n_cpu=8
)

# Run LDA topic modeling
# Test multiple topic numbers; select by model metrics
models = run_cgs_models(cistopic_obj, n_topics=[10, 20, 30, 40, 50],
                        n_cpu=8, n_iter=300, random_state=42)

# Select best model (highest coherence, lowest perplexity)
from pycisTopic.lda_models import evaluate_models
model = evaluate_models(models, return_model=True)
cistopic_obj.add_LDA_model(model)
```

## SCENIC+ Workflow

**Goal:** Assemble enhancer-driven gene regulatory networks (eRegulons) linking transcription factors to their target enhancers and downstream genes from paired scRNA+scATAC data.

**Approach:** Create a SCENICPLUS object from preprocessed scRNA-seq AnnData and cisTopic ATAC object, then run the complete pipeline which performs motif enrichment, region-to-gene linking, and eRegulon assembly.

```python
import scenicplus
from scenicplus.scenicplus_class import SCENICPLUS
from scenicplus.wrappers.run_scenicplus import run_scenicplus

# Prepare SCENIC+ object
scplus_obj = SCENICPLUS(
    adata_rna,           # scRNA-seq AnnData
    cistopic_obj,        # cisTopic object with topics
    menr=None            # will be computed
)

# Run complete SCENIC+ pipeline
# This performs: motif enrichment, region-to-gene linking, eRegulon assembly
run_scenicplus(
    scplus_obj,
    variable=['GeneExpressionLevel'],
    species='hsapiens',
    assembly='hg38',
    tf_file='allTFs_hg38.txt',
    save_path='scenicplus_output/',
    biomart_host='http://www.ensembl.org',
    upstream=[1000, 150000],    # search space upstream of TSS
    downstream=[1000, 150000],  # search space downstream of TSS
    calculate_TF_eGRN_correlation=True,
    calculate_DEGs_DARs=True,
    export_to_loom_file=True,
    export_to_UCSC_file=True,
    n_cpu=8
)
```

## eRegulon Interpretation

**Goal:** Summarize the discovered eRegulons to identify the most active transcription factors and distinguish activating from repressive regulation.

**Approach:** Parse the eRegulon metadata table to count regions and target genes per TF, then separate activating (+) and repressive (-) eRegulons by their name suffix.

```python
import pandas as pd

# eRegulons link TFs -> enhancers -> target genes
eregulons = scplus_obj.uns['eRegulon_metadata']
print(f'Found {eregulons["Region_signature_name"].nunique()} eRegulons')

# Summarize eRegulons
ereg_summary = eregulons.groupby('TF').agg(
    n_regions=('Region', 'nunique'),
    n_genes=('Gene', 'nunique')
).sort_values('n_genes', ascending=False)
print(ereg_summary.head(20))

# Activating vs repressive eRegulons
# (+) suffix = activating (TF expression positively correlated with targets)
# (-) suffix = repressive (TF expression negatively correlated with targets)
activating = eregulons[eregulons['Region_signature_name'].str.contains(r'\(\+\)')]
repressive = eregulons[eregulons['Region_signature_name'].str.contains(r'\(-\)')]
print(f'Activating: {activating["TF"].nunique()}, Repressive: {repressive["TF"].nunique()}')
```

### eRegulon Activity Scoring

```python
from scenicplus.eregulon_enrichment import score_eRegulons

# Score eRegulon activity per cell (like AUCell but for eRegulons)
score_eRegulons(scplus_obj, ranking_key='eRegulon_AUC', enrichment_type='region')
score_eRegulons(scplus_obj, ranking_key='eRegulon_AUC_gene', enrichment_type='gene')

# eRegulon AUC matrix
ereg_auc = scplus_obj.uns['eRegulon_AUC']['Gene_based'].copy()
```

## Visualization

```python
import matplotlib.pyplot as plt
from scenicplus.plotting.dotplot import heatmap_dotplot

# eRegulon activity heatmap by cell type
heatmap_dotplot(
    scplus_obj,
    size_matrix=scplus_obj.uns['eRegulon_AUC']['Gene_based'],
    color_matrix=scplus_obj.uns['eRegulon_AUC']['Region_based'],
    group_variable='cell_type',
    save='eregulon_dotplot.pdf'
)

# Network visualization of top eRegulon
from scenicplus.plotting.grn_plot import plot_eRegulon

plot_eRegulon(scplus_obj, tf_name='PAX5', save='pax5_eregulon.pdf')
```

## FigR Alternative

FigR infers TF-gene regulatory interactions from paired scRNA+scATAC without topic modeling.

```r
library(FigR)
library(Seurat)
library(Signac)

# Load Seurat object with RNA and ATAC assays (from Signac)
seurat_obj <- readRDS('multiome_seurat.rds')

# Run FigR
# Step 1: compute peak-gene correlations
peak_gene_cors <- runGenePeakcorr(
    ATAC.se = seurat_obj@assays$ATAC,
    RNAmat = seurat_obj@assays$RNA@data,
    genome = 'hg38',
    nCores = 8,
    p.cut = 0.05  # significance threshold for peak-gene links
)

# Step 2: infer TF-gene regulation scores (DORC-based)
fig_results <- runFigR(
    ATAC.se = seurat_obj@assays$ATAC,
    dorcTab = peak_gene_cors,
    genome = 'hg38',
    dorcMat = getDORCScores(seurat_obj@assays$ATAC, peak_gene_cors),
    rnaMat = seurat_obj@assays$RNA@data,
    nCores = 8
)

# Top regulatory TF-gene links
top_links <- fig_results[order(abs(fig_results$Score), decreasing = TRUE), ]
head(top_links, 20)
```

## Resource Requirements

| Dataset Size | RAM | Time | Notes |
|-------------|-----|------|-------|
| 5K cells | 16 GB | 2-4 hours | Feasible on laptop |
| 20K cells | 64 GB | 8-12 hours | Standard server |
| 50K+ cells | 128 GB+ | 24+ hours | Subsample or use HPC |

## Related Skills

- scenic-regulons - RNA-only regulon inference with pySCENIC
- single-cell/scatac-analysis - scATAC-seq preprocessing with Signac/ArchR
- atac-seq/atac-peak-calling - Peak calling for region universe
- chip-seq/motif-analysis - Motif enrichment for TF identification

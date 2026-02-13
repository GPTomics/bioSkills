---
name: bio-gene-regulatory-networks-scenic-regulons
description: Infer gene regulatory networks and identify transcription factor regulons from single-cell RNA-seq data using pySCENIC. Discovers co-expression modules with GRNBoost2, prunes by cis-regulatory motif enrichment, and scores regulon activity per cell with AUCell. Use when identifying transcription factor regulons, scoring TF activity in single cells, or finding master regulators of cell identity.
tool_type: python
primary_tool: pySCENIC
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, numpy 1.26+, pandas 2.2+, scanpy 1.10+, seaborn 0.13+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# SCENIC Regulons

**"Identify transcription factor regulons from my scRNA-seq data"** → Run the pySCENIC three-step pipeline: infer co-expression modules with GRNBoost2, prune by cis-regulatory motif enrichment with cisTarget, and score regulon activity per cell with AUCell.
- CLI: `pyscenic grn` → `pyscenic ctx` → `pyscenic aucell`
- Python: `arboreto_with_multiprocessing.py` for GRN step (workaround for dask>=2.0)

Infer transcription factor regulons from single-cell RNA-seq with the pySCENIC three-step pipeline: GRN inference, motif enrichment, and regulon activity scoring.

## Pipeline Overview

| Step | Tool | Description |
|------|------|-------------|
| 1. GRN inference | GRNBoost2 | Co-expression modules between TFs and targets |
| 2. Regulon pruning | cisTarget | Filter by cis-regulatory motif enrichment |
| 3. Activity scoring | AUCell | Score regulon activity per cell |

## Known Issues

### Arboreto / Dask Compatibility

Native Arboreto (GRNBoost2 backend) is broken with dask >= 2.0. Use the `arboreto_with_multiprocessing.py` script bundled with pySCENIC instead. This is the recommended approach for Step 1.

### Python Version

pySCENIC is tested on Python 3.10. Create a dedicated conda environment to avoid dependency conflicts:

```bash
conda create -n scenic python=3.10
conda activate scenic
pip install pyscenic loompy
```

## Required Databases

Download ranking databases and motif annotations from the cisTarget resources page (https://resources.aertslab.org/cistarget/):

```bash
# Human hg38 ranking databases (large files, ~1.5 GB each)
# mc9nr = motif collection v9, nr = non-redundant
wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather

# Motif-to-TF annotations
wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
```

## Step 1: GRN Inference with GRNBoost2

```python
import os
import glob
import pickle
import pandas as pd
import numpy as np
from pyscenic.utils import load_tf_names
from arboreto.utils import load_tf_names as arb_load_tf_names

# Load expression data (loom format is standard for SCENIC)
import loompy
ds = loompy.connect('filtered.loom')
expr_matrix = pd.DataFrame(ds[:, :], index=ds.ra.Gene, columns=ds.ca.CellID).T
ds.close()

# Load TF list (human or mouse)
tf_names = load_tf_names('allTFs_hg38.txt')
```

### Using arboreto_with_multiprocessing.py (Recommended)

```bash
# Run from command line -- avoids dask compatibility issues entirely
python arboreto_with_multiprocessing.py \
    filtered.loom \
    allTFs_hg38.txt \
    --method grnboost2 \
    --output adj.tsv \
    --num_workers 8 \
    --seed 42
```

### Python API (if dask < 2.0)

```python
from arboreto.algo import grnboost2

adjacencies = grnboost2(expr_matrix, tf_names=tf_names, verbose=True)
adjacencies.to_csv('adj.tsv', sep='\t', index=False)
```

## Step 2: Regulon Pruning with cisTarget

**Goal:** Filter the raw co-expression modules to retain only TF-target links supported by cis-regulatory motif enrichment near target gene promoters.

**Approach:** Load cisTarget ranking databases and motif annotations, run prune2df to test each TF's targets for upstream motif enrichment, then convert the pruned results into regulon objects.

```python
from pyscenic.prune import prune2df, df2regulons
from ctxcore.rnkdb import FeatherRankingDatabase

# Load ranking databases
db_fnames = glob.glob('*.genes_vs_motifs.rankings.feather')
dbs = [FeatherRankingDatabase(fname) for fname in db_fnames]

# Load motif annotations
motif_annotations_fname = 'motifs-v9-nr.hgnc-m0.001-o0.0.tbl'

adjacencies = pd.read_csv('adj.tsv', sep='\t')

# Prune: only keep TF-target links supported by cis-regulatory motifs
df = prune2df(dbs, adjacencies, motif_annotations_fname)

regulons = df2regulons(df)

with open('regulons.pkl', 'wb') as f:
    pickle.dump(regulons, f)

print(f'Found {len(regulons)} regulons')
for reg in sorted(regulons, key=lambda r: -len(r))[:10]:
    print(f'  {reg.name}: {len(reg)} targets')
```

### CLI Alternative (Steps 1-2)

```bash
# Step 1: GRN inference
pyscenic grn filtered.loom allTFs_hg38.txt -o adj.tsv --num_workers 8

# Step 2: cisTarget pruning
pyscenic ctx adj.tsv \
    hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather \
    --annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
    --expression_mtx_fname filtered.loom \
    --output reg.csv \
    --num_workers 8
```

## Step 3: AUCell Activity Scoring

**Goal:** Score the activity of each regulon in every individual cell to create a cell-by-regulon activity matrix for downstream analysis.

**Approach:** Rank genes by expression within each cell, then use AUCell to compute the area under the recovery curve for each regulon's gene set, producing an AUC score that reflects regulon activity independent of expression magnitude.

```python
from pyscenic.aucell import aucell
import loompy

ds = loompy.connect('filtered.loom')
expr_matrix = pd.DataFrame(ds[:, :], index=ds.ra.Gene, columns=ds.ca.CellID).T
ds.close()

with open('regulons.pkl', 'rb') as f:
    regulons = pickle.load(f)

# Score regulon activity per cell using AUCell
# auc_threshold: fraction of ranked genes to consider (default 0.05 = top 5%)
auc_mtx = aucell(expr_matrix, regulons, auc_threshold=0.05, num_workers=8)

auc_mtx.to_csv('auc_matrix.csv')
print(f'Scored {auc_mtx.shape[1]} regulons across {auc_mtx.shape[0]} cells')
```

### CLI Alternative

```bash
pyscenic aucell filtered.loom reg.csv \
    --output scenic_output.loom \
    --num_workers 8
```

## Interpreting Results

### Regulon Specificity Score (RSS)

```python
from pyscenic.rss import regulon_specificity_scores

# RSS identifies regulons enriched in specific cell types
# Requires cell type labels
cell_types = pd.read_csv('cell_types.csv', index_col=0)['cell_type']
rss = regulon_specificity_scores(auc_mtx, cell_types)

# Top regulons per cell type
for ct in rss.columns:
    top_regs = rss[ct].sort_values(ascending=False).head(5)
    print(f'\n{ct}:')
    for reg, score in top_regs.items():
        print(f'  {reg}: {score:.3f}')
```

### Binary Regulon Activity

```python
from pyscenic.binarization import binarize

# Binarize AUC scores (on/off per cell)
# Uses bimodal distribution fitting to set thresholds
binary_mtx, thresholds = binarize(auc_mtx)

# Fraction of cells with active regulon per cluster
cluster_activity = binary_mtx.groupby(cell_types).mean()
```

## Visualization

```python
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

adata = sc.read_h5ad('clustered.h5ad')
adata.obsm['X_aucell'] = auc_mtx.loc[adata.obs_names].values

# Regulon activity on UMAP
sc.pl.umap(adata, color=['CEBPB(+)', 'SPI1(+)', 'PAX5(+)'], cmap='viridis')

# Heatmap of top regulons per cell type
top_regulons = rss.apply(lambda x: x.nlargest(3).index.tolist()).explode().unique()
sns.clustermap(auc_mtx[top_regulons].groupby(cell_types).mean().T,
               cmap='viridis', figsize=(10, 8), z_score=0)
plt.savefig('regulon_heatmap.pdf', bbox_inches='tight')
```

## Performance Tips

| Tip | Details |
|-----|---------|
| Subsample for GRN | Use 5000-10000 cells for Step 1; regulons transfer to full dataset |
| Use CLI for Step 1 | `arboreto_with_multiprocessing.py` avoids dask issues |
| Parallelize | All three steps accept `--num_workers` |
| Prefilter genes | Remove genes expressed in < 3 cells or < 1% of cells |
| Loom format | Standard input format; convert from h5ad with `loompy` |

## Related Skills

- multiomics-grn - Enhancer-driven GRNs from paired scRNA+scATAC
- coexpression-networks - Bulk co-expression network analysis with WGCNA
- single-cell/clustering - Cluster cells before regulon analysis
- single-cell/preprocessing - QC and normalization of scRNA-seq data

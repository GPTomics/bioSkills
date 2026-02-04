---
name: bio-single-cell-splicing
description: Analyzes alternative splicing at single-cell resolution using BRIE2 for probabilistic PSI estimation or leafcutter2 for cluster-based analysis with NMD detection. Identifies cell-type-specific splicing patterns. Use when analyzing isoform usage in scRNA-seq or finding splicing differences between cell populations.
tool_type: python
primary_tool: BRIE2
---

# Single-Cell Splicing Analysis

Analyze alternative splicing at single-cell resolution.

## Tool Selection

| Tool | Approach | Strengths |
|------|----------|-----------|
| BRIE2 | Probabilistic PSI | Handles sparsity, regulatory features |
| leafcutter2 | Intron clustering | NMD detection, novel junctions |

Note: Avoid Whippet.jl (Julia 1.6.7 only, incompatible with Julia 1.9+)

## BRIE2 Analysis

```python
import brie
import scanpy as sc
import anndata as ad

# Load single-cell data
adata = sc.read_h5ad('scrnaseq.h5ad')

# Prepare splicing events from annotation
# BRIE2 uses pre-defined splicing events
brie.preprocessing.get_events(
    gtf_file='annotation.gtf',
    out_file='splicing_events.gff3'
)

# Count reads for splicing events from BAM files
# Requires cell barcodes and UMIs
brie.preprocessing.count(
    bam_file='possorted_genome_bam.bam',
    gff_file='splicing_events.gff3',
    out_dir='brie_counts/',
    cell_file='barcodes.tsv'  # Filtered cell barcodes
)

# Load BRIE count data
adata_splice = brie.read_h5ad('brie_counts/brie_count.h5ad')

# Run BRIE2 model for PSI estimation
# Uses variational inference for probabilistic estimates
brie.fit(
    adata_splice,
    layer='raw',
    n_epochs=400,
    batch_size=512
)

# PSI estimates stored in adata_splice.layers['Psi']
# Uncertainty in adata_splice.layers['Psi_var']
```

## Cell-Type Specific Splicing

```python
import numpy as np
import pandas as pd

# Add cell type annotations
adata_splice.obs['cell_type'] = adata.obs['cell_type']

# Calculate mean PSI per cell type
cell_types = adata_splice.obs['cell_type'].unique()
psi_matrix = adata_splice.layers['Psi']

mean_psi = pd.DataFrame(index=adata_splice.var_names)
for ct in cell_types:
    mask = adata_splice.obs['cell_type'] == ct
    mean_psi[ct] = np.nanmean(psi_matrix[mask, :], axis=0)

# Find cell-type specific splicing events
# Events with high variance across cell types
psi_var = mean_psi.var(axis=1)
variable_events = psi_var.nlargest(100)
print('Top variable splicing events:')
print(variable_events)
```

## leafcutter2 Analysis

```python
import subprocess

# leafcutter2 (April 2025): Adds NMD-inducing splicing detection

# Step 1: Extract junctions from BAM
# Works with 10X BAMs with cell barcodes
subprocess.run([
    'python', 'scripts/bam2junc.py',
    '-b', 'possorted_genome_bam.bam',
    '-o', 'junctions/',
    '--cb_tag', 'CB',  # Cell barcode tag
    '--umi_tag', 'UB'   # UMI tag
], check=True)

# Step 2: Cluster introns
subprocess.run([
    'python', 'clustering/leafcutter_cluster.py',
    '-j', 'junction_files.txt',
    '-o', 'leafcutter_sc',
    '-m', '10',  # Min reads per junction
    '-l', '500000'  # Max intron length
], check=True)

# Step 3: Differential splicing between clusters
# Pseudobulk approach for statistical power
subprocess.run([
    'Rscript', 'scripts/leafcutter_ds.R',
    'leafcutter_sc_perind_numers.counts.gz',
    'groups.txt',
    '-o', 'differential_splicing',
    '-e', 'annotation_exons.txt.gz'
], check=True)
```

## Pseudobulk Approach

```python
import pandas as pd
import numpy as np

# For better statistical power, aggregate cells by type
def pseudobulk_junctions(junction_counts, cell_metadata, groupby='cell_type'):
    '''Aggregate junction counts by cell group.'''
    groups = cell_metadata.groupby(groupby).groups

    pseudobulk = {}
    for group, cells in groups.items():
        cell_mask = junction_counts.index.isin(cells)
        pseudobulk[group] = junction_counts.loc[cell_mask].sum()

    return pd.DataFrame(pseudobulk)

# Run differential splicing on pseudobulk
# Use leafcutter or rMATS on aggregated counts
```

## Interpretation Considerations

| Challenge | Mitigation |
|-----------|------------|
| Sparse data | BRIE2 probabilistic model, pseudobulk |
| Low reads per cell | Aggregate similar cells |
| 3' bias (10X) | Use 5' kit or full-length methods |
| Doublets | Filter before splicing analysis |

## Quality Thresholds

| Metric | Recommendation |
|--------|----------------|
| Min cells per event | >= 50 with reads |
| Min reads per junction | >= 5 per cell with coverage |
| PSI confidence | Variance < 0.1 |

## Related Skills

- single-cell/preprocessing - QC before splicing analysis
- single-cell/clustering - Cell type annotation
- splicing-quantification - Bulk RNA-seq comparison

---
name: bio-single-cell-perturb-seq
description: Analyze Perturb-seq and CROP-seq CRISPR screening data integrated with scRNA-seq. Use when identifying gene function through pooled genetic perturbations in single cells.
tool_type: python
primary_tool: Pertpy
---

# Perturb-seq Analysis

## Load and Annotate Perturbations

```python
import scanpy as sc
import pertpy as pt

adata = sc.read_h5ad('perturb_seq.h5ad')

# Assign guide identities
# Requires guide barcode in obs or separate file
adata.obs['perturbation'] = guide_assignments
```

## Pertpy Analysis

```python
# Differential expression per perturbation
de = pt.tl.PseudobulkDE(adata)
de.fit(groupby='perturbation', control='non-targeting')
results = de.results()

# Perturbation signatures
ps = pt.tl.PerturbationSignature(adata)
ps.compute(groupby='perturbation')
```

## Mixscape (Seurat)

```r
library(Seurat)
library(SeuratObject)

# Classify perturbed vs non-perturbed cells
seurat <- RunMixscape(
    seurat,
    assay = 'RNA',
    labels = 'gene',
    nt.class.name = 'NT'
)

# Get perturbation scores
seurat <- CalcPerturbSig(
    seurat,
    assay = 'RNA',
    nt.ctrl = 'NT'
)
```

## Guide Assignment

```python
# From guide capture library
import pandas as pd

guides = pd.read_csv('guide_calls.csv')
adata.obs = adata.obs.merge(
    guides[['cell_barcode', 'guide_id', 'target_gene']],
    left_index=True, right_on='cell_barcode'
)
```

## Related Skills

- **single-cell/basic-analysis** - scRNA-seq preprocessing
- **single-cell/differential-expression** - DE concepts
- **single-cell/integration** - Multi-sample integration

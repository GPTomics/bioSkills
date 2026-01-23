---
name: bio-single-cell-lineage-tracing
description: Reconstruct cell lineage trees from CRISPR barcode tracing or mitochondrial mutations. Use when studying clonal dynamics, cell fate decisions, or developmental trajectories.
tool_type: python
primary_tool: Cassiopeia
---

# Lineage Tracing Analysis

## Cassiopeia Tree Reconstruction

```python
import cassiopeia as cas

# Load character matrix (barcode x cell)
# Values: mutation states at each barcode site
tree = cas.data.CassiopeiaTree(
    character_matrix=char_matrix,
    cell_meta=cell_metadata
)

# Reconstruct tree
solver = cas.solver.VanillaGreedySolver()
solver.solve(tree)
```

## From CRISPR Barcodes

```python
# Parse barcode sequences
barcodes = cas.pp.call_alleles(
    alignment_file='aligned_barcodes.bam',
    reference='barcode_reference.fa'
)

# Build character matrix
char_matrix = cas.pp.convert_alleles_to_character_matrix(barcodes)
```

## CoSpar for Clonal Dynamics

```python
import cospar as cs

adata = cs.read_h5ad('lineage_traced.h5ad')

# Infer transition probabilities
cs.tl.infer_Tmap(adata, smooth_array=[15, 10, 5])

# Fate probabilities
cs.tl.fate_map(adata, source='HSC', sink='Monocyte')
```

## Mitochondrial Lineage (MitoTracing)

```python
# Use mtDNA mutations as natural barcodes
import mito_lineage

variants = mito_lineage.call_variants(adata, min_cells=10)
tree = mito_lineage.build_tree(variants)
```

## Related Skills

- **single-cell/trajectory-analysis** - Pseudotime inference
- **single-cell/basic-analysis** - Preprocessing
- **phylogenetics** - Tree concepts

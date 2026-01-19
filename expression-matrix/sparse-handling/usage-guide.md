# Sparse Matrix Handling Usage Guide

## Overview

Sparse matrices store only non-zero values, making them memory-efficient for count data where many genes have zero counts across samples. This is essential for single-cell data and large bulk RNA-seq datasets.

## When to Use Sparse

| Data Type | Typical Sparsity | Use Sparse? |
|-----------|------------------|-------------|
| Bulk RNA-seq | 10-30% zeros | Sometimes |
| Single-cell | 70-95% zeros | Yes |
| Proteomics | 30-50% zeros | Maybe |
| 10X Genomics | 90%+ zeros | Yes |

Rule: Use sparse if >50% zeros and matrix is large (>10,000 genes).

## Sparse Matrix Types

```python
import scipy.sparse as sp

# CSR - Compressed Sparse Row
# Best for: row slicing, arithmetic operations
csr = sp.csr_matrix(data)

# CSC - Compressed Sparse Column
# Best for: column slicing
csc = sp.csc_matrix(data)

# COO - Coordinate format
# Best for: constructing sparse matrices, format conversion
coo = sp.coo_matrix(data)

# lil - List of Lists
# Best for: incrementally building matrix
lil = sp.lil_matrix((nrows, ncols))
```

## Complete Workflow

```python
import pandas as pd
import numpy as np
import scipy.sparse as sp
import anndata as ad

# 1. Load data
counts = pd.read_csv('counts.tsv', sep='\t', index_col=0)

# 2. Check sparsity
sparsity = (counts == 0).sum().sum() / counts.size
print(f'Sparsity: {sparsity:.1%}')

# 3. Convert to sparse if beneficial
if sparsity > 0.5:
    sparse_counts = sp.csr_matrix(counts.values)
    print(f'Memory: {counts.values.nbytes/1e6:.1f}MB -> {sparse_counts.data.nbytes/1e6:.1f}MB')
else:
    print('Matrix not sparse enough for conversion')

# 4. Create AnnData for downstream analysis
adata = ad.AnnData(
    X=sparse_counts if sparsity > 0.5 else counts.values,
    obs=pd.DataFrame(index=counts.columns),
    var=pd.DataFrame(index=counts.index)
)

# 5. Save
adata.write_h5ad('counts.h5ad', compression='gzip')
```

## Common Operations

### Filtering
```python
# Filter genes with low counts
row_sums = np.array(sparse_matrix.sum(axis=1)).flatten()
keep_genes = row_sums >= 10
filtered = sparse_matrix[keep_genes, :]

# Filter samples with low counts
col_sums = np.array(sparse_matrix.sum(axis=0)).flatten()
keep_samples = col_sums >= 1000
filtered = sparse_matrix[:, keep_samples]
```

### Normalization
```python
# CPM normalization
def sparse_cpm(X):
    lib_sizes = np.array(X.sum(axis=0)).flatten()
    return X.multiply(1e6 / lib_sizes)

# Log transform
def sparse_log1p(X):
    X = X.copy()
    X.data = np.log1p(X.data)
    return X

normalized = sparse_log1p(sparse_cpm(sparse_matrix))
```

### Statistics
```python
# Mean per gene
gene_means = np.array(sparse_matrix.mean(axis=1)).flatten()

# Variance per gene (requires dense for each row)
def sparse_var(X):
    mean = np.array(X.mean(axis=1)).flatten()
    sq_mean = np.array(X.multiply(X).mean(axis=1)).flatten()
    return sq_mean - mean**2

gene_vars = sparse_var(sparse_matrix)

# Number of non-zeros per gene
nnz_per_gene = np.array((sparse_matrix != 0).sum(axis=1)).flatten()
```

## Integration with Analysis Tools

### Scanpy
```python
import scanpy as sc

# Scanpy handles sparse automatically
adata = sc.read_h5ad('single_cell.h5ad')
print(f'Data type: {type(adata.X)}')  # Usually scipy.sparse

# Operations maintain sparsity
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
```

### Converting for DESeq2/edgeR
```python
# These tools need dense integer counts
dense_counts = sparse_matrix.toarray().astype(int)
counts_df = pd.DataFrame(dense_counts, index=gene_names, columns=sample_names)
counts_df.to_csv('counts_for_deseq.tsv', sep='\t')
```

## Memory Management

```python
import gc

# Process in chunks for very large matrices
def process_in_chunks(sparse_matrix, chunk_size=1000, func=None):
    results = []
    for i in range(0, sparse_matrix.shape[0], chunk_size):
        chunk = sparse_matrix[i:i+chunk_size, :]
        if func:
            chunk = func(chunk)
        results.append(chunk)
        gc.collect()  # Free memory
    return sp.vstack(results)

# Monitor memory
import psutil
print(f'Memory usage: {psutil.Process().memory_info().rss / 1e9:.1f} GB')
```

## Troubleshooting

### Out of Memory When Converting
```python
# Don't convert to dense
# Bad:
dense = sparse_matrix.toarray()  # May fail on large matrices

# Good: work with sparse directly
result = sparse_matrix.sum(axis=0)  # Stays sparse
```

### Slow Operations
```python
# Use appropriate format
# For row operations: CSR
# For column operations: CSC
sparse_csc = sparse_matrix.tocsc()  # Convert for column operations
col_slice = sparse_csc[:, 0:10]  # Fast column slice
```

### Indexing Returns Matrix
```python
# Single element indexing returns a matrix, not scalar
val = sparse_matrix[0, 0]  # Returns matrix
val = sparse_matrix[0, 0].item()  # Returns scalar

# Or for multiple elements
vals = np.array(sparse_matrix[0, :].todense()).flatten()
```

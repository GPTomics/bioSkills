# Sparse Matrix Handling - Usage Guide

## Overview
Work with sparse matrices (CSR, CSC, COO) for memory-efficient storage and operations on count data with many zero values, especially single-cell RNA-seq data.

## Prerequisites
```bash
pip install scipy pandas numpy anndata scanpy
```

## Quick Start
Tell your AI agent what you want to do:
- "Convert my count matrix to sparse format to save memory"
- "Calculate gene means and variances on my sparse matrix"
- "Normalize my sparse single-cell counts"

## Example Prompts
### Conversion
> "Check the sparsity of my count matrix and convert to sparse if beneficial"

> "Convert my dense matrix to CSR format for row operations"

### Operations
> "Filter genes with fewer than 10 counts across all samples from my sparse matrix"

> "Calculate CPM normalization on my sparse matrix without converting to dense"

### Integration
> "Create an AnnData object with my sparse counts and save as h5ad"

> "Convert my sparse matrix to dense for DESeq2 export"

## What the Agent Will Do
1. Assess matrix sparsity and memory footprint
2. Choose appropriate sparse format (CSR for rows, CSC for columns)
3. Convert to sparse if sparsity >50% and matrix is large
4. Perform operations while maintaining sparse format where possible
5. Convert back to dense only when required by downstream tools

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

csr = sp.csr_matrix(data)   # Best for: row slicing, arithmetic operations
csc = sp.csc_matrix(data)   # Best for: column slicing
coo = sp.coo_matrix(data)   # Best for: constructing sparse matrices, format conversion
lil = sp.lil_matrix((nrows, ncols))  # Best for: incrementally building matrix
```

## Complete Workflow

```python
import pandas as pd
import numpy as np
import scipy.sparse as sp
import anndata as ad

counts = pd.read_csv('counts.tsv', sep='\t', index_col=0)

sparsity = (counts == 0).sum().sum() / counts.size
print(f'Sparsity: {sparsity:.1%}')

if sparsity > 0.5:
    sparse_counts = sp.csr_matrix(counts.values)
    print(f'Memory: {counts.values.nbytes/1e6:.1f}MB -> {sparse_counts.data.nbytes/1e6:.1f}MB')
else:
    print('Matrix not sparse enough for conversion')

adata = ad.AnnData(
    X=sparse_counts if sparsity > 0.5 else counts.values,
    obs=pd.DataFrame(index=counts.columns),
    var=pd.DataFrame(index=counts.index)
)

adata.write_h5ad('counts.h5ad', compression='gzip')
```

## Common Operations

### Filtering
```python
row_sums = np.array(sparse_matrix.sum(axis=1)).flatten()
keep_genes = row_sums >= 10
filtered = sparse_matrix[keep_genes, :]

col_sums = np.array(sparse_matrix.sum(axis=0)).flatten()
keep_samples = col_sums >= 1000
filtered = sparse_matrix[:, keep_samples]
```

### Normalization
```python
def sparse_cpm(X):
    lib_sizes = np.array(X.sum(axis=0)).flatten()
    return X.multiply(1e6 / lib_sizes)

def sparse_log1p(X):
    X = X.copy()
    X.data = np.log1p(X.data)
    return X

normalized = sparse_log1p(sparse_cpm(sparse_matrix))
```

### Statistics
```python
gene_means = np.array(sparse_matrix.mean(axis=1)).flatten()

def sparse_var(X):
    mean = np.array(X.mean(axis=1)).flatten()
    sq_mean = np.array(X.multiply(X).mean(axis=1)).flatten()
    return sq_mean - mean**2

gene_vars = sparse_var(sparse_matrix)

nnz_per_gene = np.array((sparse_matrix != 0).sum(axis=1)).flatten()
```

## Integration with Analysis Tools

### Scanpy
```python
import scanpy as sc

adata = sc.read_h5ad('single_cell.h5ad')
print(f'Data type: {type(adata.X)}')

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
```

### Converting for DESeq2/edgeR
```python
dense_counts = sparse_matrix.toarray().astype(int)
counts_df = pd.DataFrame(dense_counts, index=gene_names, columns=sample_names)
counts_df.to_csv('counts_for_deseq.tsv', sep='\t')
```

## Memory Management

```python
import gc

def process_in_chunks(sparse_matrix, chunk_size=1000, func=None):
    results = []
    for i in range(0, sparse_matrix.shape[0], chunk_size):
        chunk = sparse_matrix[i:i+chunk_size, :]
        if func:
            chunk = func(chunk)
        results.append(chunk)
        gc.collect()
    return sp.vstack(results)

import psutil
print(f'Memory usage: {psutil.Process().memory_info().rss / 1e9:.1f} GB')
```

## Troubleshooting

### Out of Memory When Converting
```python
dense = sparse_matrix.toarray()  # May fail on large matrices

result = sparse_matrix.sum(axis=0)  # Stays sparse - preferred
```

### Slow Operations
```python
sparse_csc = sparse_matrix.tocsc()
col_slice = sparse_csc[:, 0:10]  # Fast column slice with CSC format
```

### Indexing Returns Matrix
```python
val = sparse_matrix[0, 0]         # Returns matrix
val = sparse_matrix[0, 0].item()  # Returns scalar

vals = np.array(sparse_matrix[0, :].todense()).flatten()
```

## Tips
- Use CSR format for row operations and CSC for column operations
- Avoid converting to dense unless absolutely necessary (e.g., for DESeq2)
- Process in chunks for very large matrices to manage memory
- Scanpy and AnnData handle sparse matrices natively - let them do the work
- Monitor memory usage with psutil when working with large datasets

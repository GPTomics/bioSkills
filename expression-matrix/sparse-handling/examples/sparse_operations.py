#!/usr/bin/env python3
'''Handle sparse count matrices for single-cell data'''

import scipy.sparse as sp
import anndata as ad
import numpy as np
import pandas as pd

def dense_to_sparse(dense_matrix):
    '''Convert dense matrix to sparse CSR format'''
    return sp.csr_matrix(dense_matrix)

def sparse_to_dense(sparse_matrix):
    '''Convert sparse matrix to dense (use with caution for large matrices)'''
    return sparse_matrix.toarray()

def load_10x_mtx(path):
    '''Load 10X Genomics matrix.mtx format'''
    adata = ad.read_10x_mtx(path)
    return adata

def save_sparse_h5ad(adata, path):
    '''Save AnnData with sparse matrix to h5ad'''
    if not sp.issparse(adata.X):
        adata.X = sp.csr_matrix(adata.X)
    adata.write(path)

def sparse_row_sums(sparse_matrix):
    '''Efficient row sums for sparse matrix'''
    return np.array(sparse_matrix.sum(axis=1)).flatten()

def sparse_col_sums(sparse_matrix):
    '''Efficient column sums for sparse matrix'''
    return np.array(sparse_matrix.sum(axis=0)).flatten()

def filter_sparse_genes(adata, min_cells=3):
    '''Filter genes by minimum cell count'''
    gene_counts = sparse_col_sums(adata.X > 0)
    keep = gene_counts >= min_cells
    return adata[:, keep].copy()

def subsample_sparse(adata, n_cells=1000):
    '''Subsample cells from sparse AnnData'''
    if adata.n_obs <= n_cells:
        return adata
    idx = np.random.choice(adata.n_obs, n_cells, replace=False)
    return adata[idx].copy()

if __name__ == '__main__':
    # Example with random sparse data
    dense = np.random.poisson(0.1, (1000, 500))
    sparse = dense_to_sparse(dense)
    print(f'Dense size: {dense.nbytes / 1e6:.1f} MB')
    print(f'Sparse size: {(sparse.data.nbytes + sparse.indices.nbytes + sparse.indptr.nbytes) / 1e6:.1f} MB')
    print(f'Sparsity: {1 - sparse.nnz / np.prod(sparse.shape):.2%}')

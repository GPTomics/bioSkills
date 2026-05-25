---
name: bio-expression-matrix-sparse-handling
description: Stores and operates on sparse expression matrices for single-cell and large bulk RNA-seq, covering dgCMatrix/dgRMatrix/dgTMatrix when-each-is-fast, the dgCMatrix (CSC, R) <-> CSR (Python) implicit transpose, AnnData (cells-rows) <-> SingleCellExperiment (cells-cols) orientation flip, HDF5/h5ad vs Zarr cloud-native shift, HDF5SummarizedExperiment + DelayedArray for out-of-memory bulk, scanpy backed mode for large h5ad, the ~10-15% density crossover where dense beats sparse, 10X format proliferation (MTX vs CellRanger H5 vs h5ad), the dense-conversion memory blow-up, and Dask + Zarr for consortium-scale matrices. Use when choosing sparse format, working with single-cell-sized matrices, importing/exporting 10X, debugging R/Python interop transposes, processing matrices too large for RAM, or building cloud-native pipelines.
tool_type: python
primary_tool: scipy.sparse
---

## Version Compatibility

Reference examples tested with: numpy 1.26+, scipy 1.12+, pandas 2.2+, anndata 0.10+, scanpy 1.10+, Matrix R package 1.6+, HDF5Array 1.30+ (Bioconductor), DelayedArray 0.28+, zellkonverter 1.12+, zarr-python 2.18+, dask 2024.1+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Sparse Matrix Handling

**"Store / compute on a single-cell matrix without blowing up memory"** -> Pick the sparse format that matches the access pattern (CSC for column ops, CSR for row ops), respect the R/Python convention difference (Bioconductor stores cells in columns; AnnData stores cells in rows), and use HDF5/Zarr backed mode for matrices too large for RAM.

## The Single Most Important Modern Insight -- The R/Python interop transpose is silent and catastrophic

R Bioconductor (Seurat, SingleCellExperiment) stores cells in COLUMNS and uses dgCMatrix (Column-Compressed Sparse). Python scverse (AnnData, scanpy) stores cells in ROWS and defaults to CSR (Compressed Sparse Row). Round-tripping with `anndata2ri`, `zellkonverter`, or `rpy2` triggers a TRANSPOSE under the hood -- once per direction. Two flips silently cancel. A debugging session that converts back and forth multiple times can end up with mysteriously transposed data and no error.

| Conversion | Implicit transpose |
|------------|-------------------|
| Python CSR -> R dgCMatrix | YES (rows <-> cols) |
| AnnData `.X` (cells x genes) -> SingleCellExperiment counts (genes x cells) | YES |
| Seurat `@assays$RNA@counts` (genes x cells) -> AnnData `.X` (cells x genes) | YES |
| `scipy.sparse.csr_matrix(dense)` | None (just format conversion) |
| `csr.tocsc()` / `csc.tocsr()` | None (semantically same matrix, different layout) |

The safe pattern for one-shot conversions: file-based intermediate. `adata.write('file.h5ad')` then `zellkonverter::readH5AD('file.h5ad')` -- avoids the in-memory rpy2/reticulate gymnastics, and the file roundtrip makes orientation explicit.

For a 1M-cell single-cell matrix, the implicit transpose is non-trivial -- minutes of wall time and a temporary memory peak roughly equal to nnz x 12 bytes (CSC) or x 16 bytes (CSR with int64). Avoid unnecessary transposes by aligning the format to the consumer.

Two adjacent insights:

1. **Sparse becomes inefficient above ~10-15% density.** Sparse iteration has cache-unfriendly indirection; dense iteration is sequential. Above ~10-15% nonzero, dense is often faster for most operations even though it uses more memory.
2. **AnnData backed mode quietly differs from in-memory in important ways.** `sc.read_h5ad('file.h5ad', backed='r')` returns an AnnData where `.X` is a wrapped HDF5 dataset, read-only. Many scanpy functions silently load to memory; some functions error or hang on backed mode.

## Algorithmic Taxonomy

| Format | Layout | Fast for | Slow for |
|--------|--------|----------|----------|
| dgCMatrix (CSC, R) | i (row indices), p (col pointers), x (values) | Column slicing; per-cell ops in single-cell (cells in cols); matrix-vector with col vector | Row slicing |
| dgRMatrix (CSR, R) | j (col indices), p (row pointers), x (values) | Row slicing; per-gene ops when genes in rows | Column slicing |
| dgTMatrix (COO, R triplet) | i, j, x | Random insertion when building; reading MTX | Most operations -- convert to dgCMatrix after build |
| scipy `csc_matrix` | Same as dgCMatrix | Column ops in Python | Row ops |
| scipy `csr_matrix` | Same as dgRMatrix | Row ops (NumPy convention); matmul; sklearn defaults | Column ops |
| scipy `coo_matrix` | Triplet (i, j, data) | Construction; MTX I/O | Most ops -- convert after build |
| HDF5 (h5ad, h5) | Single-file binary chunked | Random access via chunks; compression; widely supported | Cloud / parallel writes |
| Zarr | Chunked array, per-chunk file (or S3 object) | Cloud-native; parallel writes; Dask integration | Single-file simplicity |
| DelayedArray + HDF5Array | Bioc lazy evaluation over HDF5 | Out-of-memory bulk ops in R | Speed of in-memory |

## Decision Tree by Scenario

| Scenario | Recommended approach |
|----------|---------------------|
| Single-cell (>10k cells), per-cell ops | dgCMatrix in R (Bioconductor); CSR `adata.X` in Python (scanpy) |
| Bulk RNA-seq (60-80% density) | Dense -- sparse overhead exceeds benefit |
| Single-cell pseudobulk (after donor aggregation) | Dense -- now 60-80% density typically |
| 1M+ cells, can't fit in RAM | scanpy backed mode OR HDF5SummarizedExperiment + DelayedArray |
| TCGA + GTEx + recount3 scale (100k+ samples) | HDF5Array / Zarr + Dask |
| 10X CellRanger 3.0+ output | `sc.read_10x_h5()` or `Read10X_h5()` -- the .h5 is faster than the .mtx triplet |
| Cloud-native (anndata on S3, dask compute) | Zarr |
| Local workstation, single-machine | HDF5 -- faster, more widely supported |
| Building sparse matrix incrementally | COO (dgTMatrix / coo_matrix); convert to CSC/CSR after |
| R <-> Python conversion | File-based intermediate (`adata.write` then `zellkonverter::readH5AD`); aware of the transpose |

## Check Sparsity

**Goal:** Decide whether sparse is the right format given the actual data density.

**Approach:** Compute nonzero fraction; rule of thumb is sparse > ~85% (single-cell). Below that, dense often wins.

```python
import numpy as np
import scipy.sparse as sp

def sparsity(m):
    if sp.issparse(m):
        return 1 - m.nnz / (m.shape[0] * m.shape[1])
    return (m == 0).mean()

s = sparsity(adata.X)
print(f'{s:.1%} sparse')
```

Memory math:

| Data | Format | Bytes |
|------|--------|-------|
| 30k x 100k single-cell matrix, 5% density | dgCMatrix | (5% * 3e9) * 12 bytes ~= 1.8 GB |
| Same | Dense double | 24 GB |
| 60k x 100k bulk, 70% density | dgCMatrix | 50 GB (worse than dense!) |
| Same | Dense double | 48 GB |

For single-cell (typically 90-95% sparse), sparse is essential. For bulk RNA-seq (typically 60-80% density), dense is faster and not appreciably larger.

## dgCMatrix / scipy CSC / CSR

**Goal:** Construct, query, and convert sparse matrices in the format matching the consumer's expected layout.

**Approach:** `Matrix::sparseMatrix(i, j, x, dims=...)` (R) or `scipy.sparse.csr_matrix((data, (i, j)))` (Python); preserve row/column names; convert layout (CSC <-> CSR) without changing semantics.

```python
import scipy.sparse as sp
import pandas as pd

dense_df = pd.read_csv('counts.csv', index_col=0)
sparse_csr = sp.csr_matrix(dense_df.values)
sparse_csc = sp.csc_matrix(dense_df.values)

gene_names = dense_df.index.tolist()
sample_names = dense_df.columns.tolist()

sparse_csr.tocsc()
sparse_csc.tocsr()
```

```r
library(Matrix)

dense_mat <- as.matrix(read.csv('counts.csv', row.names = 1))
sparse_dgc <- as(dense_mat, 'CsparseMatrix')

class(sparse_dgc)
rownames(sparse_dgc) <- rownames(dense_mat)
colnames(sparse_dgc) <- colnames(dense_mat)
```

dgTMatrix is best for building matrices incrementally (reading MTX, parsing per-row); convert to dgCMatrix for downstream ops:

```r
mat_t <- as(triplet_data, 'TsparseMatrix')
mat_c <- as(mat_t, 'CsparseMatrix')
```

## HDF5 vs Zarr -- The Cloud-Native Shift

HDF5 (Hierarchical Data Format 5): hierarchical, single-file binary. Random access via chunks; supports compression (gzip, blosc, lz4). On-disk format for AnnData `.h5ad`, MuData `.h5mu`, 10x Genomics `.h5`, HDF5SummarizedExperiment.

Zarr: cloud-native, chunked array storage. Each chunk is a separate file (or S3 object). Parallel-write friendly; splittable by Dask. Format used by recent AnnData (`anndata.write_zarr`), SpatialData, and large-cohort consortia.

| Criterion | HDF5 | Zarr |
|-----------|------|------|
| File structure | Single binary file | Directory of chunk files |
| Parallel writes | Limited (process-level locks) | Native |
| S3 / cloud object storage | Workarounds (h5cloud); often slow | Native; first-class |
| Compression options | gzip, blosc, lz4, szip | gzip, blosc, lz4, zstd, custom |
| Local workstation speed | Faster | Slightly slower (many small files) |
| Wide ecosystem support | Yes (mature) | Growing; modern scverse |

For local workstation work, HDF5 is faster and more widely supported. For cloud-mounted analysis (anndata on S3 with dask-distributed compute), Zarr wins because of object-storage friendliness.

## HDF5SummarizedExperiment + DelayedArray (Bioconductor)

**Goal:** Work with bulk SummarizedExperiment objects too large to fit in RAM by keeping the matrix on disk.

**Approach:** `HDF5Array` wraps an HDF5 dataset as a `DelayedArray`. Subsetting builds a delayed operation tree -- no I/O until realization. `DelayedMatrixStats` provides delayed-friendly stat functions.

```r
library(HDF5Array)
library(SummarizedExperiment)
library(DelayedMatrixStats)

se <- loadHDF5SummarizedExperiment('saved_se_dir')

s <- se[1:1000, 1:50]
row_means <- rowMeans2(assay(se))

saveHDF5SummarizedExperiment(se, 'saved_se_dir', replace = TRUE)
```

The `assay(se)` returns a DelayedMatrix backed by HDF5. DESeq2, edgeR, limma have varying levels of DelayedArray support; consult package docs before assuming all ops work in delayed mode.

For TCGA + GTEx + recount3 scale (100k+ samples, 60k genes), a dense matrix is ~48 GB (double); dgCMatrix at 70% density is ~50 GB (sparse loses). HDF5Array + chunk-aware ops keeps memory at whatever-fits-in-RAM.

## scanpy Backed Mode

**Goal:** Work with h5ad files too large for memory by loading only accessed slices on demand.

**Approach:** `sc.read_h5ad(..., backed='r')` returns an AnnData with `.X` as a wrapped HDF5 dataset. Subset operations are lazy; `.to_memory()` realizes.

```python
import scanpy as sc

adata = sc.read_h5ad('large_dataset.h5ad', backed='r')
print(f'Shape: {adata.shape}, X type: {type(adata.X)}')

t_cells = adata[adata.obs['cell_type'] == 'T_cell', :].to_memory()
```

Limitations:
- `.X` is read-only in `backed='r'`. Use `backed='r+'` for in-place updates, but only `.X` updates supported.
- `.obs` and `.var` are fully loaded -- only `.X` supports backed access.
- Very large sparse h5ad (>35 GB) can still cause memory issues even in backed mode (anndata library overhead).
- Many scanpy functions internally load to memory; check `?function` docs for backed compatibility.
- Functions like `sc.tl.pca`, `sc.pp.neighbors` typically require in-memory; subset first with `.to_memory()`.

For datasets too large for backed mode, process in chunks:

```python
import anndata as ad

def process_in_chunks(h5ad_path, chunk_size=10000, func=None):
    adata = sc.read_h5ad(h5ad_path, backed='r')
    n_cells = adata.shape[0]
    results = []
    for start in range(0, n_cells, chunk_size):
        end = min(start + chunk_size, n_cells)
        chunk = adata[start:end].to_memory()
        if func:
            chunk = func(chunk)
        results.append(chunk)
    return ad.concat(results)
```

## 10X Genomics Format Proliferation

| Format | Files | Notes |
|--------|-------|-------|
| MTX (pre-CellRanger 3.0) | `matrix.mtx` + `barcodes.tsv` + `features.tsv` (or `genes.tsv`) | Triplet format; slow to read for large matrices |
| H5 (CellRanger 3.0+) | `filtered_feature_bc_matrix.h5` | HDF5 with `/matrix/data`, `/matrix/indices`, `/matrix/indptr`, `/matrix/shape`; single file, fast |
| H5AD | `data.h5ad` | AnnData; convert on import |
| kallisto|bustools output | `output.bus` + barcode and gene mappings | `BUSpaRse` / `kb-python` |

```python
import scanpy as sc

adata = sc.read_10x_h5('filtered_feature_bc_matrix.h5')
adata = sc.read_10x_mtx('filtered_feature_bc_matrix/')
```

```r
library(Seurat)
mat <- Read10X_h5('filtered_feature_bc_matrix.h5')
mat <- Read10X(data.dir = 'filtered_feature_bc_matrix/')

library(DropletUtils)
sce <- read10xCounts('filtered_feature_bc_matrix/')
```

For 10X output, prefer the `.h5` over the MTX triplet -- typically 5-10x faster for large matrices.

## Dense Conversion -- The Memory Blow-Up

`adata.X.toarray()` (Python) or `as.matrix(seurat_obj@assays$RNA@counts)` (R) on a 30k x 100k single-cell matrix instantiates a ~24 GB dense double array. Common triggers:

- Passing sparse to a function that internally calls `as.matrix()` (older R `cor()` implementations).
- Heatmap functions (`pheatmap`, `ComplexHeatmap`) that require dense.
- ML libraries with no sparse support (some sklearn models; XGBoost requires specific sparse API).
- Plotting functions (`plot()`, `ggplot2`) called on the full matrix.

Defensive pattern: subset to a manageable gene/cell set BEFORE dense conversion.

For per-cell PCA-style ops, use sparse-aware solvers:

```python
from scipy.sparse.linalg import svds
U, s, Vt = svds(adata.X, k=50)
```

```r
library(irlba)
svd_res <- irlba(sparse_mat, nv = 50)
```

`irlba` (R) and `scipy.sparse.linalg.svds` (Python) compute truncated SVD without densifying.

## SCE vs AnnData vs MuData -- Where Bulk Fits

| Container | Library | Cells/samples | Multi-modal | Bulk fit |
|-----------|---------|---------------|-------------|----------|
| SummarizedExperiment / RangedSummarizedExperiment | Bioconductor | n/a; bulk | No | YES -- standard for DESeq2/edgeR/limma bulk |
| SingleCellExperiment (Amezquita 2020) | Bioconductor | cells in cols | Via `altExps` | scRNA-seq with spike-ins, ADT |
| AnnData | scverse/Python | cells in rows | Via layers | scRNA-seq; bulk is unusual |
| MuData | scverse/Python | cells in rows | Yes, multiple AnnData | Multi-modal scRNA + ATAC + protein |
| MultiAssayExperiment | Bioconductor | samples | Yes | R-side multi-modal analog |

Bulk RNA-seq rarely uses AnnData -- it shines on the single-cell dimensionality reduction / neighbors / clustering machinery. For bulk in R, use SummarizedExperiment. For bulk in Python, a tidy DataFrame + numpy array is usually sufficient.

## Dask + Zarr for Consortium-Scale Matrices

For TCGA + GTEx + recount3 (Wilks 2021 *Genome Biol* 22:323) or pancancer assemblies:

```python
import zarr
import dask.array as da

z = zarr.open('counts.zarr', mode='r')
da_arr = da.from_zarr(z)

col_sums = da_arr.sum(axis=0).compute()
filtered = da_arr[da_arr.sum(axis=1) > 100, :]
```

```python
import anndata as ad
adata_disk = ad.read_zarr('large_data.zarr')
```

For R: `HDF5Array` + `DelayedArray` is the equivalent, but R doesn't have a true Dask analog. BiocParallel can parallelize chunks, but lazy planning is more manual.

## Sparse Operations

```python
import numpy as np
import scipy.sparse as sp

row_sums = np.array(sparse_matrix.sum(axis=1)).flatten()
col_sums = np.array(sparse_matrix.sum(axis=0)).flatten()

keep_rows = row_sums > 10
sparse_filt = sparse_matrix[keep_rows, :]

sparse_log = sparse_matrix.copy()
sparse_log.data = np.log1p(sparse_log.data)
```

Subsetting: select genes (rows) or samples (cols) by index:

```python
gene_idx = [gene_names.index(g) for g in ['TP53', 'BRCA1', 'MYC'] if g in gene_names]
subset = sparse_matrix[gene_idx, :]
```

## CPM Normalization on Sparse

**Goal:** Apply CPM normalization without densifying.

**Approach:** Compute library sizes from column sums; broadcast scaling factors with sparse multiply for CPM; transform only the nonzero data array in-place with log1p.

```python
import numpy as np
import scipy.sparse as sp

def normalize_sparse_cpm(sparse_matrix):
    lib_sizes = np.array(sparse_matrix.sum(axis=0)).flatten()
    scaling = 1e6 / lib_sizes
    return sparse_matrix.multiply(scaling)

def log1p_inplace(sparse_matrix):
    out = sparse_matrix.copy()
    out.data = np.log1p(out.data)
    return out

cpm = normalize_sparse_cpm(adata.X)
log_cpm = log1p_inplace(cpm)
```

After log-transformation, sparsity is PRESERVED (log1p(0) = 0). After CPM with pseudocount, zeros become nonzero -- check sparsity and convert to dense if density drops below ~15%.

## Save / Load Sparse Matrices

```python
import scipy.sparse as sp
import numpy as np

sp.save_npz('counts_sparse.npz', sparse_matrix)
loaded = sp.load_npz('counts_sparse.npz')

np.savez('counts_with_meta.npz',
    data    = sparse_matrix.data,
    indices = sparse_matrix.indices,
    indptr  = sparse_matrix.indptr,
    shape   = sparse_matrix.shape,
    genes   = np.array(gene_names),
    samples = np.array(sample_names))
```

For interop and durable storage, prefer h5ad or zarr:

```python
adata.write_h5ad('counts.h5ad')
adata.write_zarr('counts.zarr')
```

## Per-Method Failure Modes

### Implicit transpose in R/Python conversion

**Trigger:** AnnData with cells in rows passed to a SingleCellExperiment workflow that expects cells in cols; downstream `colSums` returns gene-level totals.

**Mechanism:** AnnData stores cells in rows; SCE in cols. The conversion auto-transposes ONCE per direction; two roundtrips silently restore.

**Symptom:** Per-cell stats look like per-gene stats; QC plots have wrong axes.

**Fix:** Use file-based intermediate (`adata.write('file.h5ad')`; `zellkonverter::readH5AD('file.h5ad')`). Always verify dimensions and orientation after conversion.

### Dense conversion blew up memory

**Trigger:** `as.matrix(seurat_obj@assays$RNA@counts)` on a 100k-cell dataset; R session crashes with OOM.

**Mechanism:** 100k cells x 30k genes = 3e9 entries; double precision = 24 GB.

**Symptom:** R session killed; "cannot allocate vector of size N GB".

**Fix:** Don't densify the full matrix. Subset to genes/cells of interest first. For dimensionality reduction, use `irlba::irlba()` (sparse SVD).

### scanpy backed mode silently loaded to memory

**Trigger:** `adata = sc.read_h5ad(path, backed='r')` then `sc.tl.pca(adata)`; memory spikes to dense-equivalent.

**Mechanism:** Many scanpy functions internally call `.to_memory()` because they cannot operate on backed mode. `sc.tl.pca`, `sc.pp.neighbors`, `sc.tl.umap` all materialize.

**Symptom:** OOM despite backed mode.

**Fix:** Subset first (`adata[mask].to_memory()`), then operate. Or use a streaming-aware alternative (Dask + Zarr).

### Sparse stored where dense would be faster

**Trigger:** Bulk RNA-seq with 70% density stored as dgCMatrix; per-gene `rowVars` is slow.

**Mechanism:** Sparse iteration has cache-unfriendly indirection; above ~10-15% density, dense wins.

**Symptom:** Operations notably slower than expected; profiler shows time in sparse indexing.

**Fix:** Convert to dense for the hot path: `as.matrix(sparse_mat)` (R) or `sparse_matrix.toarray()` (Python). Memory may go up but speed improves substantially.

### 10X MTX read is slow

**Trigger:** Reading a 100k-cell 10X dataset via the MTX three-file format; takes 10+ minutes.

**Mechanism:** MTX is a text format; parsing is slow for large matrices.

**Symptom:** Long load times; user kills the process before completion.

**Fix:** Use the CellRanger H5 (`.h5`) instead -- typically 5-10x faster.

## Common errors

| Error / symptom | Cause | Fix |
|-----------------|-------|-----|
| `cannot allocate vector of size N GB` | Implicit dense conversion | Subset first; use sparse-aware solver (irlba) |
| Sparse-dense arithmetic returns `numpy.matrix` | Deprecated NumPy type from sparse+dense | `np.asarray(sparse + dense)` to force ndarray |
| `KeyError: '_index'` reading h5ad | anndata version mismatch | Update anndata; or `sc.read_h5ad(..., backed=None)` |
| Empty rows/cols after sparse subset | Subset removed all data | Verify the index list; cross-check sample/gene names |
| Backed mode AnnData crash on `sc.tl.umap` | Function not backed-compatible | `.to_memory()` on the subset first |
| Per-cell totals look wrong after R<->Python conversion | Implicit transpose | Verify dimensions; use file-based intermediate |
| CSR (Python) <-> dgCMatrix (R) treated as same | Convention difference | They're transposes of each other; verify shape and a known cell-gene pair |

## References

- Amezquita RA, Lun ATL, Becht E et al. 2020. Orchestrating single-cell analysis with Bioconductor. *Nat Methods* 17:137-145. doi:10.1038/s41592-019-0654-x
- Wolf FA, Angerer P, Theis FJ. 2018. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol* 19:15. doi:10.1186/s13059-017-1382-0
- Wilks C et al. 2021. recount3: summaries and queries for large-scale RNA-seq expression and splicing. *Genome Biol* 22:323. doi:10.1186/s13059-021-02533-6
- Bates D, Maechler M. 2023. Matrix: Sparse and Dense Matrix Classes and Methods. R package version 1.6-x.
- Pages H et al. 2020. HDF5Array: HDF5 backend for DelayedArray objects. Bioconductor package.
- Miles A et al. 2020. zarr-python. Python package documentation.
- Rocklin M. 2015. Dask: Parallel Computation with Blocked algorithms and Task Scheduling. *Proc Python Sci Conf.* (canonical Dask reference)
- Lachmann A et al. 2018. Massive mining of publicly available RNA-seq data from human and mouse. *Nat Commun* 9:1366. doi:10.1038/s41467-018-03751-6

## Related Skills

- counts-ingest - Reading 10X formats; building sparse matrices from quantification output
- gene-id-mapping - Var (gene) metadata in AnnData
- metadata-joins - Obs (sample) metadata in AnnData
- normalization - log1p and CPM patterns on sparse
- differential-expression/deseq2-basics - Pseudobulk aggregation makes dense
- single-cell/data-io - Single-cell file format ecosystem
- single-cell/preprocessing - Standard single-cell sparse pipeline

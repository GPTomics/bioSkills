---
name: bio-expression-matrix-counts-ingest
description: Load gene expression count matrices from various formats including CSV, TSV, featureCounts, Salmon, kallisto, and 10X. Use when importing quantification results for downstream analysis.
tool_type: python
primary_tool: pandas
---

## Version Compatibility

Reference examples tested with: pandas 2.2+, tximport 1.30+, tximeta 1.20+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Count Matrix Ingestion

## Expression Units Decision Table

| Unit | Source | Use for DE | Use for visualization | Cross-sample comparison |
|------|--------|------------|----------------------|------------------------|
| Raw counts | featureCounts, HTSeq, STAR | Yes (DE tools normalize internally) | No | No |
| Estimated counts | Salmon, kallisto (via tximport) | Yes (with offset correction) | No | No |
| TPM | Salmon, kallisto | No | Within-sample gene comparison | Partially (same composition caveat) |
| FPKM/RPKM | Cufflinks, legacy tools | No | Within-sample only | No (composition bias) |
| Normalized counts | DESeq2, edgeR | Via DE tool | Prefer VST/rlog | Yes |

For differential expression, always provide raw integer counts. DE tools (DESeq2, edgeR, limma-voom) model count distributions and handle normalization internally.

## Transcript-Level vs Gene-Level Quantification

Salmon and kallisto quantify at the transcript level. Naive summation of transcript counts to gene level introduces **length bias**: if treatment causes a shift from short to long isoforms, the gene appears upregulated even if the number of mRNA molecules is unchanged. Use tximport in R (or equivalent offset handling in Python) to correct for this.

For detailed tximport workflows, see rna-quantification/tximport-workflow. Key decisions:
- `countsFromAbundance='no'` (default): raw estimated counts with length offsets for DESeq2 (recommended)
- `countsFromAbundance='lengthScaledTPM'`: required for limma-voom, which cannot use offset matrices
- `countsFromAbundance='scaledTPM'`: prevents differential transcript usage from creating spurious DGE

## Basic CSV/TSV Loading

**Goal:** Load a gene expression count matrix from a delimited text file into a pandas DataFrame.

**Approach:** Read CSV/TSV with gene IDs as the row index, handling comment lines if present.

**"Load my count matrix"** → Read a delimited file into a DataFrame with genes as rows and samples as columns.

```python
import pandas as pd

# TSV with gene IDs as first column
counts = pd.read_csv('counts.tsv', sep='\t', index_col=0)

# CSV with header
counts = pd.read_csv('counts.csv', index_col=0)

# Skip comment lines
counts = pd.read_csv('counts.txt', sep='\t', index_col=0, comment='#')
```

## featureCounts Output

**Goal:** Parse featureCounts output into a clean count matrix by stripping metadata columns.

**Approach:** Skip the 6 annotation columns (Chr, Start, End, Strand, Length) and clean BAM path suffixes from column names.

```python
import pandas as pd

# featureCounts format has 6 metadata columns before counts
fc = pd.read_csv('featurecounts.txt', sep='\t', comment='#')
counts = fc.set_index('Geneid').iloc[:, 5:]  # Skip Chr, Start, End, Strand, Length
counts.columns = [c.replace('.bam', '').split('/')[-1] for c in counts.columns]
```

## Salmon Quant Files

**Goal:** Combine per-sample Salmon quantification files into a single count or TPM matrix.

**Approach:** Iterate over quant directories, extract the desired column from each quant.sf, and merge into a DataFrame.

```python
import pandas as pd
from pathlib import Path

def load_salmon_quants(quant_dirs, column='NumReads'):
    '''Load multiple Salmon quant.sf files into a count matrix.'''
    dfs = {}
    for qdir in quant_dirs:
        sample = Path(qdir).name
        sf = pd.read_csv(f'{qdir}/quant.sf', sep='\t', index_col=0)
        dfs[sample] = sf[column]
    return pd.DataFrame(dfs)

# Usage
quant_dirs = ['salmon_out/sample1', 'salmon_out/sample2', 'salmon_out/sample3']
counts = load_salmon_quants(quant_dirs, column='NumReads')
tpm = load_salmon_quants(quant_dirs, column='TPM')
```

## STAR ReadsPerGene Files

**Goal:** Load STAR gene-level count output, selecting the correct strandedness column.

**Approach:** STAR outputs 4 columns per gene: gene ID, unstranded counts, sense-strand counts, antisense-strand counts. Selecting the wrong strand column silently halves the counts.

```python
import pandas as pd
from pathlib import Path

def load_star_genecounts(filepaths, strandedness='reverse'):
    '''Load STAR ReadsPerGene.out.tab files.

    strandedness: 'unstranded' (col 1), 'forward' (col 2), 'reverse' (col 3)
    For Illumina TruSeq stranded libraries, use 'reverse'.
    '''
    col_map = {'unstranded': 1, 'forward': 2, 'reverse': 3}
    col_idx = col_map[strandedness]
    dfs = {}
    for fp in filepaths:
        sample = Path(fp).name.replace('_ReadsPerGene.out.tab', '')
        df = pd.read_csv(fp, sep='\t', header=None, index_col=0)
        dfs[sample] = df.iloc[4:, col_idx - 1]  # Skip first 4 summary rows
    return pd.DataFrame(dfs)
```

Strandedness verification: compare total counts in columns 2 vs 3 across samples. The column with much larger totals corresponds to the correct strand. If roughly equal, the library is unstranded (use column 1).

## HTSeq Count Files

**Goal:** Load per-sample HTSeq count output into a combined matrix.

**Approach:** Read each file (two columns: gene ID and count), skipping the 5 summary lines at the end (prefixed with `__`), then merge.

```python
import pandas as pd
from pathlib import Path

def load_htseq_counts(filepaths):
    '''Load multiple HTSeq count files into a matrix.'''
    dfs = {}
    for fp in filepaths:
        sample = Path(fp).stem.replace('_counts', '')
        df = pd.read_csv(fp, sep='\t', header=None, index_col=0, names=['gene', 'count'])
        df = df[~df.index.str.startswith('__')]  # Remove __no_feature, __ambiguous, etc.
        dfs[sample] = df['count']
    return pd.DataFrame(dfs)
```

## kallisto Abundance Files

**Goal:** Combine per-sample kallisto abundance files into a single count or TPM matrix.

**Approach:** Read each abundance.tsv, extract the target column, and assemble into a DataFrame keyed by sample name.

```python
import pandas as pd
from pathlib import Path

def load_kallisto_quants(abundance_files, column='est_counts'):
    '''Load multiple kallisto abundance.tsv files.'''
    dfs = {}
    for f in abundance_files:
        sample = Path(f).parent.name
        ab = pd.read_csv(f, sep='\t', index_col=0)
        dfs[sample] = ab[column]
    return pd.DataFrame(dfs)

# Usage
files = ['kallisto_out/sample1/abundance.tsv', 'kallisto_out/sample2/abundance.tsv']
counts = load_kallisto_quants(files, column='est_counts')
tpm = load_kallisto_quants(files, column='tpm')
```

## 10X Genomics Sparse Matrix

**Goal:** Load single-cell count data from 10X Genomics format (MTX directory or H5 file).

**Approach:** Use scanpy to read the sparse matrix with gene and barcode metadata into an AnnData object.

```python
import scanpy as sc

# Load 10X directory (contains matrix.mtx, genes.tsv/features.tsv, barcodes.tsv)
adata = sc.read_10x_mtx('filtered_feature_bc_matrix/')

# Load 10X H5 file
adata = sc.read_10x_h5('filtered_feature_bc_matrix.h5')

# Convert to dense DataFrame if needed
counts = adata.to_df()
```

## AnnData H5AD Files

**Goal:** Load preprocessed expression data stored in the AnnData H5AD format.

**Approach:** Read the H5AD file with scanpy and access the count matrix, raw layer, or metadata as needed.

```python
import anndata as ad
import scanpy as sc

# Load h5ad
adata = sc.read_h5ad('data.h5ad')

# Access count matrix
counts = adata.to_df()  # Dense DataFrame
sparse_counts = adata.X  # Sparse matrix (if stored sparse)

# Access raw counts if normalized data is in .X
raw_counts = adata.raw.to_adata().to_df()
```

## RDS Files (from R)

**Goal:** Load R-serialized count data into Python.

**Approach:** Use pyreadr to read RDS files directly into a pandas DataFrame.

```python
import pyreadr

# Read RDS file
result = pyreadr.read_r('counts.rds')
counts = result[None]  # Access the data

# For Seurat objects, use anndata2ri or convert in R first
```

## Combine Multiple Files

**Goal:** Merge per-sample count files into a single genes-by-samples matrix.

**Approach:** Glob matching files, read the first data column from each, and concatenate into a DataFrame.

```python
import pandas as pd
from pathlib import Path

def combine_count_files(file_pattern, index_col=0, sep='\t'):
    '''Combine multiple count files into one matrix.'''
    files = sorted(Path('.').glob(file_pattern))
    dfs = {}
    for f in files:
        sample = f.stem.replace('_counts', '')
        dfs[sample] = pd.read_csv(f, sep=sep, index_col=index_col).iloc[:, 0]
    return pd.DataFrame(dfs)

# Usage
counts = combine_count_files('counts/*_counts.tsv')
```

## Filter Low-Count Genes

**Goal:** Remove genes with insufficient expression to reduce noise and multiple testing burden.

**Approach:** Apply expression thresholds appropriate for the downstream DE tool. edgeR requires explicit filtering; DESeq2 handles it internally but benefits from a speed-only pre-filter.

### edgeR: Design-Aware Filtering (Recommended)

```r
library(edgeR)

# filterByExpr uses smallest group size from design matrix
# Requires CPM above threshold in at least n samples (n = smallest group)
keep <- filterByExpr(y, design=model.matrix(~condition, data=metadata))
y <- y[keep, , keep.lib.sizes=FALSE]
```

Forgetting `filterByExpr` in edgeR inflates the multiple testing burden because edgeR's quasi-likelihood pipeline does not have DESeq2's automatic independent filtering.

### DESeq2: Speed Pre-Filter

```r
# Manual pre-filter for speed only -- independent filtering at results() handles statistics
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
```

### Python

```python
# Group-aware filter: require min_counts in at least min_samples
min_counts, min_samples = 10, 3
expressed = (counts >= min_counts).sum(axis=1) >= min_samples
counts_filtered = counts.loc[expressed]

# For min_samples, use the smallest experimental group size
# e.g., if 3 controls and 5 treated, use min_samples=3
```

## Handle Gene ID Versions

**Goal:** Strip Ensembl version suffixes for consistent ID matching across datasets.

**Approach:** Split on the dot separator and keep only the stable identifier prefix.

```python
# Remove Ensembl version numbers (ENSG00000123456.12 -> ENSG00000123456)
counts.index = counts.index.str.split('.').str[0]

# Or keep as-is for compatibility
```

## Save Count Matrix

**Goal:** Export a count matrix in formats suitable for downstream tools or archival.

**Approach:** Write as TSV, compressed TSV, or AnnData H5AD depending on the target workflow.

```python
# Save as TSV
counts.to_csv('count_matrix.tsv', sep='\t')

# Save as compressed
counts.to_csv('count_matrix.tsv.gz', sep='\t', compression='gzip')

# Save as AnnData
import anndata as ad
adata = ad.AnnData(counts)
adata.write_h5ad('counts.h5ad')
```

## R Loading Equivalents

**Goal:** Load count data in R from common formats including featureCounts and Salmon/kallisto via tximport.

**Approach:** Use base R read functions for alignment-based counts, or tximport with a tx2gene mapping for pseudo-alignment outputs (Salmon/kallisto). tximport handles length-bias correction automatically.

```r
# Basic CSV/TSV
counts <- read.csv('counts.csv', row.names=1)
counts <- read.delim('counts.tsv', row.names=1)

# featureCounts
fc <- read.delim('featurecounts.txt', comment.char='#', row.names=1)
counts <- fc[, 6:ncol(fc)]

# tximport for Salmon/kallisto (RECOMMENDED over manual loading)
# Use tx2gene for gene-level summarization with length-offset correction
library(tximport)
tx2gene <- read.csv('tx2gene.csv')  # columns: TXNAME, GENEID
files <- file.path('salmon_out', samples, 'quant.sf')
names(files) <- samples
txi <- tximport(files, type='salmon', tx2gene=tx2gene)

# For DESeq2: use DESeqDataSetFromTximport (handles offsets natively)
# For limma-voom: use countsFromAbundance='lengthScaledTPM' in tximport
# See rna-quantification/tximport-workflow for full details
```

## Related Skills

- rna-quantification/featurecounts-counting - Generate featureCounts output
- rna-quantification/alignment-free-quant - Generate Salmon/kallisto output
- rna-quantification/tximport-workflow - Import Salmon/kallisto with length-offset correction
- expression-matrix/normalization - Normalize counts for downstream analysis
- expression-matrix/sparse-handling - Memory-efficient storage
- expression-matrix/gene-id-mapping - Convert gene identifiers

# Count Matrix Ingestion - Usage Guide

## Overview
Load gene expression count matrices from various quantification tools (featureCounts, Salmon, kallisto, STAR, HTSeq, 10X) into pandas DataFrames for downstream analysis.

## Prerequisites
```bash
pip install pandas numpy scipy
```

## Quick Start
Tell your AI agent what you want to do:
- "Load my featureCounts output file into a pandas DataFrame"
- "Combine Salmon quantification files from multiple samples"
- "Import STAR gene counts with reverse strand specificity"

## Example Prompts
### Loading Specific Formats
> "Load the featureCounts file at counts.txt and clean up the sample names"

> "Read all Salmon quant.sf files from the salmon_quants directory into a combined count matrix"

> "Import STAR ReadsPerGene.out.tab files with reverse strandedness"

### Quality Checks
> "Load my count matrix and show me basic QC statistics"

> "Check my count matrix for duplicate gene IDs and sum them"

### Format Detection
> "Figure out what delimiter my counts file uses and load it"

## What the Agent Will Do
1. Identify the input format based on file structure or user specification
2. Load the data with appropriate parsing (skip headers, set index, clean sample names)
3. Run basic quality checks (shape, library sizes, zero counts, NaN values)
4. Optionally convert to sparse format for memory efficiency if matrix is highly sparse

## Common Formats

| Source | Format | Key Columns |
|--------|--------|-------------|
| featureCounts | TSV | Geneid + 5 meta cols + counts |
| Salmon | quant.sf per sample | NumReads, TPM |
| kallisto | abundance.tsv per sample | est_counts, tpm |
| STAR | ReadsPerGene.out.tab | gene_id, unstranded/stranded counts |
| 10X | matrix.mtx + features/barcodes | Sparse format |
| HTSeq | TSV | gene_id, count |

## Complete Loading Pipeline

```python
import pandas as pd
import numpy as np
from pathlib import Path

class CountMatrixLoader:
    @staticmethod
    def from_featurecounts(filepath):
        df = pd.read_csv(filepath, sep='\t', comment='#')
        counts = df.set_index('Geneid').iloc[:, 5:]
        counts.columns = [c.replace('.bam', '').split('/')[-1] for c in counts.columns]
        return counts

    @staticmethod
    def from_salmon_dir(base_dir):
        base = Path(base_dir)
        samples = [d.name for d in base.iterdir() if d.is_dir() and (d / 'quant.sf').exists()]
        dfs = {}
        for sample in samples:
            sf = pd.read_csv(base / sample / 'quant.sf', sep='\t', index_col=0)
            dfs[sample] = sf['NumReads']
        return pd.DataFrame(dfs)

    @staticmethod
    def from_star_genecounts(filepaths, strandedness='reverse'):
        col_map = {'unstranded': 1, 'forward': 2, 'reverse': 3}
        col_idx = col_map[strandedness]
        dfs = {}
        for fp in filepaths:
            sample = Path(fp).name.replace('_ReadsPerGene.out.tab', '')
            df = pd.read_csv(fp, sep='\t', header=None, index_col=0)
            dfs[sample] = df.iloc[4:, col_idx - 1]
        return pd.DataFrame(dfs)

counts = CountMatrixLoader.from_featurecounts('counts.txt')
counts = CountMatrixLoader.from_salmon_dir('salmon_quants/')
```

## Quality Checks After Loading

```python
def check_count_matrix(counts):
    print(f'Shape: {counts.shape[0]} genes x {counts.shape[1]} samples')
    print(f'Total counts per sample:\n{counts.sum().describe()}')
    print(f'Genes with zero counts: {(counts.sum(axis=1) == 0).sum()}')
    print(f'Any NaN values: {counts.isna().any().any()}')
    print(f'Library sizes range: {counts.sum().min():.0f} - {counts.sum().max():.0f}')
    return counts

counts = check_count_matrix(counts)
```

## Handling Different Organisms

| Organism | Ensembl Format | Example |
|----------|----------------|---------|
| Human | ENSG | ENSG00000141510 |
| Mouse | ENSMUSG | ENSMUSG00000059552 |
| Zebrafish | ENSDARG | ENSDARG00000002354 |
| Fly | FBgn | FBgn0000008 |
| Worm | WBGene | WBGene00000001 |

```python
sample_ids = counts.index[:5].tolist()
print(sample_ids)

if counts.index.str.startswith('ENSG').any():
    print('Human Ensembl gene IDs')
elif counts.index.str.startswith('ENSMUSG').any():
    print('Mouse Ensembl gene IDs')
```

## Troubleshooting

### Duplicate Gene IDs
```python
if counts.index.duplicated().any():
    print(f'Duplicate IDs: {counts.index.duplicated().sum()}')
    counts = counts.groupby(counts.index).sum()
```

### Missing Samples
```python
expected = ['sample1', 'sample2', 'sample3']
actual = counts.columns.tolist()
missing = set(expected) - set(actual)
if missing:
    print(f'Missing samples: {missing}')
```

### Wrong Delimiter
```python
import csv
with open('counts.txt', 'r') as f:
    dialect = csv.Sniffer().sniff(f.read(1024))
    print(f'Detected delimiter: {repr(dialect.delimiter)}')
```

## Tips
- Always check the matrix shape and library sizes after loading to catch parsing errors early
- Remove version suffixes from Ensembl IDs (e.g., ENSG00000141510.15 -> ENSG00000141510) before downstream analysis
- Use sparse matrices for single-cell data or bulk RNA-seq with >90% zeros
- Batch query multiple samples at once rather than loading one at a time for better performance

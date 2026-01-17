# Count Matrix Ingestion Usage Guide

## Overview

This skill covers loading gene expression count matrices from various quantification tools and formats. Proper loading and initial processing is essential before differential expression or other downstream analyses.

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
        '''Load featureCounts output.'''
        df = pd.read_csv(filepath, sep='\t', comment='#')
        counts = df.set_index('Geneid').iloc[:, 5:]
        counts.columns = [c.replace('.bam', '').split('/')[-1] for c in counts.columns]
        return counts

    @staticmethod
    def from_salmon_dir(base_dir):
        '''Load Salmon quants from directory structure.'''
        base = Path(base_dir)
        samples = [d.name for d in base.iterdir() if d.is_dir() and (d / 'quant.sf').exists()]
        dfs = {}
        for sample in samples:
            sf = pd.read_csv(base / sample / 'quant.sf', sep='\t', index_col=0)
            dfs[sample] = sf['NumReads']
        return pd.DataFrame(dfs)

    @staticmethod
    def from_star_genecounts(filepaths, strandedness='reverse'):
        '''Load STAR ReadsPerGene.out.tab files.'''
        col_map = {'unstranded': 1, 'forward': 2, 'reverse': 3}
        col_idx = col_map[strandedness]
        dfs = {}
        for fp in filepaths:
            sample = Path(fp).name.replace('_ReadsPerGene.out.tab', '')
            df = pd.read_csv(fp, sep='\t', header=None, index_col=0)
            dfs[sample] = df.iloc[4:, col_idx - 1]  # Skip first 4 rows (summary)
        return pd.DataFrame(dfs)

# Usage
counts = CountMatrixLoader.from_featurecounts('counts.txt')
counts = CountMatrixLoader.from_salmon_dir('salmon_quants/')
```

## Handling Different Organisms

Gene IDs vary by organism and annotation source:

| Organism | Ensembl Format | Example |
|----------|----------------|---------|
| Human | ENSG | ENSG00000141510 |
| Mouse | ENSMUSG | ENSMUSG00000059552 |
| Zebrafish | ENSDARG | ENSDARG00000002354 |
| Fly | FBgn | FBgn0000008 |
| Worm | WBGene | WBGene00000001 |

```python
# Check ID format
sample_ids = counts.index[:5].tolist()
print(sample_ids)

# Common patterns
if counts.index.str.startswith('ENSG').any():
    print('Human Ensembl gene IDs')
elif counts.index.str.startswith('ENSMUSG').any():
    print('Mouse Ensembl gene IDs')
```

## Quality Checks After Loading

```python
def check_count_matrix(counts):
    '''Run basic QC on count matrix.'''
    print(f'Shape: {counts.shape[0]} genes x {counts.shape[1]} samples')
    print(f'Total counts per sample:\n{counts.sum().describe()}')
    print(f'Genes with zero counts: {(counts.sum(axis=1) == 0).sum()}')
    print(f'Any NaN values: {counts.isna().any().any()}')
    print(f'Library sizes range: {counts.sum().min():.0f} - {counts.sum().max():.0f}')
    return counts

counts = check_count_matrix(counts)
```

## Memory Optimization

For large matrices, use sparse representation:

```python
import scipy.sparse as sp

# Convert to sparse if >90% zeros
sparsity = (counts == 0).sum().sum() / counts.size
print(f'Matrix sparsity: {sparsity:.1%}')

if sparsity > 0.9:
    sparse_counts = sp.csr_matrix(counts.values)
    print(f'Dense size: {counts.values.nbytes / 1e6:.1f} MB')
    print(f'Sparse size: ~{sparse_counts.data.nbytes / 1e6:.1f} MB')
```

## Batch Loading Script

```python
#!/usr/bin/env python
import pandas as pd
import argparse
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description='Load count matrices')
    parser.add_argument('--format', choices=['featurecounts', 'salmon', 'kallisto', 'tsv'])
    parser.add_argument('--input', required=True, help='Input file or directory')
    parser.add_argument('--output', required=True, help='Output file')
    args = parser.parse_args()

    if args.format == 'featurecounts':
        counts = load_featurecounts(args.input)
    elif args.format == 'salmon':
        counts = load_salmon_dir(args.input)
    elif args.format == 'kallisto':
        counts = load_kallisto_dir(args.input)
    else:
        counts = pd.read_csv(args.input, sep='\t', index_col=0)

    counts.to_csv(args.output, sep='\t')
    print(f'Saved {counts.shape[0]} genes x {counts.shape[1]} samples to {args.output}')

if __name__ == '__main__':
    main()
```

## Troubleshooting

### Duplicate Gene IDs
```python
# Check for duplicates
if counts.index.duplicated().any():
    print(f'Duplicate IDs: {counts.index.duplicated().sum()}')
    # Sum duplicates
    counts = counts.groupby(counts.index).sum()
```

### Missing Samples
```python
# Check expected vs actual samples
expected = ['sample1', 'sample2', 'sample3']
actual = counts.columns.tolist()
missing = set(expected) - set(actual)
if missing:
    print(f'Missing samples: {missing}')
```

### Wrong Delimiter
```python
# Auto-detect delimiter
import csv
with open('counts.txt', 'r') as f:
    dialect = csv.Sniffer().sniff(f.read(1024))
    print(f'Detected delimiter: {repr(dialect.delimiter)}')
```

# Metadata Joins Usage Guide

## Overview

Properly joining sample metadata with count matrices is critical for differential expression and other analyses. Mismatched samples or missing metadata can cause analysis failures or incorrect results.

## Metadata Structure

A well-structured metadata file should have:
- Sample IDs as row index (matching count matrix columns)
- Experimental factors as columns
- No missing values in critical columns

```csv
sample_id,condition,batch,sex,age
sample1,control,batch1,M,45
sample2,control,batch1,F,52
sample3,treated,batch2,M,48
sample4,treated,batch2,F,51
```

## Complete Workflow

```python
import pandas as pd
import anndata as ad

class DataPreparation:
    def __init__(self, counts_file, metadata_file):
        self.counts = pd.read_csv(counts_file, sep='\t', index_col=0)
        self.metadata = pd.read_csv(metadata_file, index_col=0)
        self.harmonize()

    def harmonize(self):
        '''Match samples between counts and metadata.'''
        count_samples = set(self.counts.columns)
        meta_samples = set(self.metadata.index)
        common = count_samples & meta_samples

        if len(common) == 0:
            raise ValueError('No matching samples between counts and metadata')

        dropped_counts = count_samples - common
        dropped_meta = meta_samples - common

        if dropped_counts:
            print(f'Dropping {len(dropped_counts)} samples not in metadata')
        if dropped_meta:
            print(f'Dropping {len(dropped_meta)} metadata rows not in counts')

        self.counts = self.counts[sorted(common)]
        self.metadata = self.metadata.loc[sorted(common)]

    def validate(self, required_columns):
        '''Check metadata completeness.'''
        for col in required_columns:
            if col not in self.metadata.columns:
                raise ValueError(f'Required column missing: {col}')
            if self.metadata[col].isna().any():
                raise ValueError(f'Missing values in column: {col}')
        return True

    def to_anndata(self):
        '''Export as AnnData.'''
        adata = ad.AnnData(X=self.counts.T)
        adata.obs = self.metadata
        return adata

    def to_deseq_files(self, prefix):
        '''Export files for DESeq2.'''
        self.counts.to_csv(f'{prefix}_counts.tsv', sep='\t')
        self.metadata.to_csv(f'{prefix}_metadata.csv')

# Usage
prep = DataPreparation('counts.tsv', 'metadata.csv')
prep.validate(['condition', 'batch'])
adata = prep.to_anndata()
```

## Handling Common Issues

### Sample Name Variations
```python
def normalize_sample_names(names):
    '''Standardize sample names.'''
    normalized = []
    for n in names:
        n = str(n)
        n = n.replace('.bam', '')
        n = n.replace('_sorted', '')
        n = n.split('/')[-1]  # Remove path
        normalized.append(n)
    return normalized

counts.columns = normalize_sample_names(counts.columns)
metadata.index = normalize_sample_names(metadata.index)
```

### Categorical Variables
```python
# Ensure categorical columns are properly typed
categorical_cols = ['condition', 'batch', 'sex']
for col in categorical_cols:
    if col in metadata.columns:
        metadata[col] = pd.Categorical(metadata[col])

# Set reference level for differential expression
metadata['condition'] = pd.Categorical(
    metadata['condition'],
    categories=['control', 'treated'],
    ordered=True
)
```

### Missing Values
```python
# Check for missing values
print(metadata.isna().sum())

# Options for handling:
# 1. Drop samples with any missing
metadata_clean = metadata.dropna()

# 2. Fill with defaults
metadata['batch'] = metadata['batch'].fillna('unknown')

# 3. Impute numeric columns
metadata['age'] = metadata['age'].fillna(metadata['age'].median())
```

## R Workflow

```r
library(DESeq2)

# Load and validate
counts <- read.delim('counts.tsv', row.names=1, check.names=FALSE)
metadata <- read.csv('metadata.csv', row.names=1)

# Match samples
stopifnot(all(colnames(counts) %in% rownames(metadata)))
metadata <- metadata[colnames(counts), , drop=FALSE]

# Set factor levels (control as reference)
metadata$condition <- factor(metadata$condition, levels=c('control', 'treated'))

# Validate
stopifnot(all(colnames(counts) == rownames(metadata)))
stopifnot(!any(is.na(metadata$condition)))

# Create DESeq object
dds <- DESeqDataSetFromMatrix(
    countData=round(as.matrix(counts)),
    colData=metadata,
    design=~batch + condition
)
```

## Multi-factor Designs

```python
# Combine factors for complex designs
metadata['group'] = metadata['condition'] + '_' + metadata['timepoint']

# Or keep separate for interaction models
# In DESeq2: design = ~condition + timepoint + condition:timepoint
```

```r
# R version
metadata$group <- paste(metadata$condition, metadata$timepoint, sep='_')

# Or use interaction in design
dds <- DESeqDataSetFromMatrix(
    countData=counts,
    colData=metadata,
    design=~condition * timepoint
)
```

## Batch Information

Always include batch information if available:

```python
# Add sequencing batch
metadata['seq_batch'] = metadata['sample_id'].apply(
    lambda x: 'batch1' if x.startswith('S1') else 'batch2'
)

# Add processing date if available
# This helps identify technical variation
```

## Quality Control Integration

```python
# Add QC metrics from alignment/counting
qc_metrics = pd.read_csv('alignment_qc.csv', index_col=0)
metadata = metadata.join(qc_metrics, how='left')

# Flag low-quality samples
metadata['pass_qc'] = (
    (metadata['uniquely_mapped_pct'] > 70) &
    (metadata['total_reads'] > 1e6) &
    (metadata['assigned_pct'] > 50)
)

# Filter to passing samples
counts = counts.loc[:, metadata[metadata['pass_qc']].index]
metadata = metadata[metadata['pass_qc']]
```

## Export for Various Tools

```python
def export_for_deseq2(counts, metadata, prefix):
    '''Export for DESeq2 analysis in R.'''
    # Counts must be integers
    counts_int = counts.round().astype(int)
    counts_int.to_csv(f'{prefix}_counts.tsv', sep='\t')
    metadata.to_csv(f'{prefix}_coldata.csv')

def export_for_edger(counts, metadata, prefix):
    '''Export for edgeR analysis in R.'''
    counts.to_csv(f'{prefix}_counts.txt', sep='\t')
    # edgeR uses group from metadata
    metadata.to_csv(f'{prefix}_design.csv')

def export_for_limma(counts, metadata, prefix):
    '''Export for limma-voom.'''
    # Similar to DESeq2
    counts.to_csv(f'{prefix}_counts.txt', sep='\t')
    metadata.to_csv(f'{prefix}_targets.csv')
```

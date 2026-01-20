# Metadata Joins - Usage Guide

## Overview
Join sample metadata with count matrices ensuring proper sample alignment, validating metadata completeness, and exporting for downstream differential expression analysis tools.

## Prerequisites
```bash
pip install pandas anndata
```

```r
install.packages("BiocManager")
BiocManager::install("DESeq2")
```

## Quick Start
Tell your AI agent what you want to do:
- "Join my metadata file with my count matrix"
- "Validate that all samples have condition and batch information"
- "Export my data for DESeq2 analysis"

## Example Prompts
### Basic Joining
> "Load my counts.tsv and metadata.csv and join them together"

> "Match samples between my count matrix and sample metadata, keeping only common samples"

### Validation
> "Check that all samples have values for condition, batch, and sex columns"

> "Find samples that are in my counts but missing from metadata"

### Export
> "Export my data as files ready for DESeq2 analysis in R"

> "Create an AnnData object with my counts and metadata"

### Sample Name Issues
> "My sample names in counts have .bam suffixes but metadata doesn't - harmonize them"

## What the Agent Will Do
1. Load count matrix and metadata files
2. Harmonize sample names (clean suffixes, match formats)
3. Identify common samples and report any dropped
4. Validate required columns exist and have no missing values
5. Set appropriate factor levels for categorical variables
6. Export in the format needed for downstream tools

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
        for col in required_columns:
            if col not in self.metadata.columns:
                raise ValueError(f'Required column missing: {col}')
            if self.metadata[col].isna().any():
                raise ValueError(f'Missing values in column: {col}')
        return True

    def to_anndata(self):
        adata = ad.AnnData(X=self.counts.T)
        adata.obs = self.metadata
        return adata

    def to_deseq_files(self, prefix):
        self.counts.to_csv(f'{prefix}_counts.tsv', sep='\t')
        self.metadata.to_csv(f'{prefix}_metadata.csv')

prep = DataPreparation('counts.tsv', 'metadata.csv')
prep.validate(['condition', 'batch'])
adata = prep.to_anndata()
```

## Handling Common Issues

### Sample Name Variations
```python
def normalize_sample_names(names):
    normalized = []
    for n in names:
        n = str(n)
        n = n.replace('.bam', '')
        n = n.replace('_sorted', '')
        n = n.split('/')[-1]
        normalized.append(n)
    return normalized

counts.columns = normalize_sample_names(counts.columns)
metadata.index = normalize_sample_names(metadata.index)
```

### Categorical Variables
```python
categorical_cols = ['condition', 'batch', 'sex']
for col in categorical_cols:
    if col in metadata.columns:
        metadata[col] = pd.Categorical(metadata[col])

metadata['condition'] = pd.Categorical(metadata['condition'], categories=['control', 'treated'], ordered=True)
```

### Missing Values
```python
print(metadata.isna().sum())

metadata_clean = metadata.dropna()
metadata['batch'] = metadata['batch'].fillna('unknown')
metadata['age'] = metadata['age'].fillna(metadata['age'].median())
```

## R Workflow

```r
library(DESeq2)

counts <- read.delim('counts.tsv', row.names=1, check.names=FALSE)
metadata <- read.csv('metadata.csv', row.names=1)

stopifnot(all(colnames(counts) %in% rownames(metadata)))
metadata <- metadata[colnames(counts), , drop=FALSE]

metadata$condition <- factor(metadata$condition, levels=c('control', 'treated'))

stopifnot(all(colnames(counts) == rownames(metadata)))
stopifnot(!any(is.na(metadata$condition)))

dds <- DESeqDataSetFromMatrix(countData=round(as.matrix(counts)), colData=metadata, design=~batch + condition)
```

## Multi-factor Designs

```python
metadata['group'] = metadata['condition'] + '_' + metadata['timepoint']
```

```r
metadata$group <- paste(metadata$condition, metadata$timepoint, sep='_')

dds <- DESeqDataSetFromMatrix(countData=counts, colData=metadata, design=~condition * timepoint)
```

## Quality Control Integration

```python
qc_metrics = pd.read_csv('alignment_qc.csv', index_col=0)
metadata = metadata.join(qc_metrics, how='left')

metadata['pass_qc'] = (
    (metadata['uniquely_mapped_pct'] > 70) &
    (metadata['total_reads'] > 1e6) &
    (metadata['assigned_pct'] > 50)
)

counts = counts.loc[:, metadata[metadata['pass_qc']].index]
metadata = metadata[metadata['pass_qc']]
```

## Export for Various Tools

```python
def export_for_deseq2(counts, metadata, prefix):
    counts_int = counts.round().astype(int)
    counts_int.to_csv(f'{prefix}_counts.tsv', sep='\t')
    metadata.to_csv(f'{prefix}_coldata.csv')

def export_for_edger(counts, metadata, prefix):
    counts.to_csv(f'{prefix}_counts.txt', sep='\t')
    metadata.to_csv(f'{prefix}_design.csv')

def export_for_limma(counts, metadata, prefix):
    counts.to_csv(f'{prefix}_counts.txt', sep='\t')
    metadata.to_csv(f'{prefix}_targets.csv')
```

## Tips
- Always verify sample order matches between counts and metadata before analysis
- Set factor reference levels explicitly (e.g., 'control' before 'treated') for interpretable results
- Include batch information whenever available to account for technical variation
- Validate metadata completeness early to avoid cryptic errors in downstream tools
- Use descriptive column names and document the meaning of categorical levels

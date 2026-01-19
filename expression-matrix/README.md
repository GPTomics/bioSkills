# expression-matrix

## Overview

Load, manipulate, and annotate gene expression count matrices. Covers reading various formats (CSV, TSV, H5AD, RDS), sparse matrix handling for large datasets, gene ID mapping between databases, and joining sample metadata.

**Tool type:** mixed | **Primary tools:** pandas, anndata, scanpy, biomaRt, pyensembl

## Skills

| Skill | Description |
|-------|-------------|
| counts-ingest | Load count matrices from CSV, TSV, featureCounts, Salmon, 10X formats |
| sparse-handling | Work with sparse matrices for memory-efficient storage |
| gene-id-mapping | Convert between Ensembl, Entrez, HGNC, and gene symbols |
| metadata-joins | Merge sample metadata with count matrices |

## Example Prompts

- "Load my featureCounts output into a dataframe"
- "Read a 10X sparse matrix"
- "Convert my count matrix to sparse format"
- "Map Ensembl IDs to gene symbols"
- "Convert Entrez IDs to Ensembl"
- "Join sample metadata with my count matrix"
- "Load an AnnData h5ad file"
- "Combine multiple count files into one matrix"
- "Filter my count matrix by expressed genes"
- "Add gene annotations to my matrix"

## Requirements

```bash
# Python
pip install pandas numpy scipy anndata scanpy pyensembl

# For gene ID mapping
pip install mygene gseapy
pyensembl install --release 110 --species human
```

```r
# R/Bioconductor
BiocManager::install(c('biomaRt', 'AnnotationDbi', 'org.Hs.eg.db', 'org.Mm.eg.db'))
```

## Related Skills

- **rna-quantification** - Generate count matrices
- **single-cell** - Single-cell specific data handling
- **differential-expression** - Downstream DE analysis
- **pathway-analysis** - Functional enrichment analysis

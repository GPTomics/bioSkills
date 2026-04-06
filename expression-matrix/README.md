# expression-matrix

## Overview

Load, normalize, manipulate, and annotate gene expression count matrices. Covers reading various formats (CSV, TSV, H5AD, RDS, 10X), normalization and transformation for different downstream tasks, sparse matrix handling for large datasets, gene ID mapping between databases, and joining sample metadata with experimental design guidance.

**Tool type:** mixed | **Primary tools:** pandas, anndata, DESeq2, edgeR, scanpy, biomaRt, pyensembl

## Skills

| Skill | Description |
|-------|-------------|
| counts-ingest | Load count matrices from CSV, TSV, featureCounts, Salmon, kallisto, STAR, HTSeq, 10X formats |
| normalization | Normalize and transform counts (TMM, RLE, VST, rlog, TPM) for DE, visualization, or clustering |
| sparse-handling | Work with sparse matrices for memory-efficient storage and backed mode for large datasets |
| gene-id-mapping | Convert between Ensembl, Entrez, HGNC, UniProt; cross-species orthologs; tx2gene construction |
| metadata-joins | Merge sample metadata with count matrices; sample swap detection; experimental design guidance |

## Example Prompts

- "Load my featureCounts output into a dataframe"
- "Import Salmon quant files using tximport with proper length-offset correction"
- "Which normalization should I use for PCA visualization vs DE analysis?"
- "Apply VST transformation to my count matrix for clustering"
- "Convert my Ensembl IDs to gene symbols"
- "Map gene IDs from human to mouse orthologs"
- "Join sample metadata with my count matrix and check for sample swaps"
- "Read a 10X sparse matrix in backed mode"
- "Filter my count matrix by expressed genes using edgeR's filterByExpr"
- "Build a tx2gene mapping from my GTF file"

## Requirements

```bash
# Python
pip install pandas numpy scipy anndata scanpy pyensembl mygene pydeseq2

# For gene ID mapping
pyensembl install --release 110 --species human
```

```r
# R/Bioconductor
BiocManager::install(c('DESeq2', 'edgeR', 'limma', 'tximport', 'tximeta',
                        'biomaRt', 'AnnotationDbi', 'org.Hs.eg.db', 'org.Mm.eg.db',
                        'scran', 'EDASeq', 'cqn'))
```

## Related Skills

- **rna-quantification** - Generate count matrices and QC
- **single-cell** - Single-cell specific data handling
- **differential-expression** - Downstream DE analysis
- **pathway-analysis** - Functional enrichment analysis

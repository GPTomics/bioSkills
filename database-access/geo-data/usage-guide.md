# GEO Data - Usage Guide

## Overview

This skill enables AI agents to help you search and access gene expression data from NCBI GEO (Gene Expression Omnibus), including microarray and RNA-seq datasets.

## Prerequisites

```bash
pip install biopython

# Optional for full data parsing
pip install GEOparse
```

## Quick Start

Tell your AI agent what you want to do:

- "Find RNA-seq datasets for breast cancer in humans"
- "Get the GEO series associated with PMID 35412348"
- "Find the SRA accessions for GSE123456"
- "Download the expression matrix for this GEO dataset"

## Example Prompts

### Find Datasets
> "Search GEO for human single-cell RNA-seq datasets from 2024"

### By Disease/Condition
> "Find GEO datasets studying Alzheimer's disease in mouse"

### From Publication
> "What GEO datasets are associated with this PubMed article?"

### Get SRA Data
> "I want to download the raw reads from GSE123456 - find the SRA accessions"

### Download Data
> "How do I download the processed expression data from GSE123456?"

## What the Agent Will Do

1. Search the GDS database using appropriate filters
2. Return dataset summaries with accessions, organisms, and sample counts
3. Link to SRA data if available
4. Provide download options (FTP, GEOparse)

## GEO Record Types

- **GSE** (Series) - Complete experiment with multiple samples
- **GSM** (Sample) - Individual sample
- **GPL** (Platform) - Array or sequencing platform
- **GDS** (DataSet) - Curated, normalized dataset

## Tips

- Specify organism when searching
- Use "RNA-seq" or "microarray" to filter technology
- For raw sequencing data, link GSE to SRA
- Series matrix files contain processed expression values
- Use GEOparse for full metadata and data parsing

# GTF/GFF Handling - Usage Guide

## Overview

GTF (Gene Transfer Format) and GFF3 (General Feature Format) are standard formats for gene annotations. This skill covers parsing, querying, and converting these files using gtfparse, gffutils, and gffread.

## When to Use This Skill

- You need to extract gene coordinates from annotation files
- You want to convert between GTF and GFF3 formats
- You need to get promoters, exons, or introns for analysis
- You want to query annotations programmatically

## Installation

```bash
# gffread (CLI - format conversion)
conda install -c bioconda gffread

# Python packages
pip install gtfparse gffutils
```

## Key Concepts

### Coordinate Systems

Both GTF and GFF3 use **1-based, inclusive** coordinates:
- First base is position 1
- Start and end are both included
- BED uses 0-based, half-open (subtract 1 from start when converting)

### GTF vs GFF3

| Use GTF when | Use GFF3 when |
|--------------|---------------|
| Working with Ensembl/GENCODE | Working with NCBI RefSeq |
| Need transcript structure | Need explicit parent-child |
| Standard RNA-seq pipelines | Complex nested features |

## Basic Workflow

### 1. Download Annotations

```bash
# GENCODE GTF (human)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
gunzip gencode.v44.annotation.gtf.gz

# Ensembl GFF3 (human)
wget https://ftp.ensembl.org/pub/release-110/gff3/homo_sapiens/Homo_sapiens.GRCh38.110.gff3.gz
gunzip Homo_sapiens.GRCh38.110.gff3.gz
```

### 2. Parse and Extract

```python
import gtfparse

df = gtfparse.read_gtf('gencode.v44.annotation.gtf')

# Get protein-coding genes
coding_genes = df[(df['feature'] == 'gene') &
                  (df['gene_type'] == 'protein_coding')]
print(f'Found {len(coding_genes)} protein-coding genes')
```

### 3. Convert to BED

```python
# Convert to BED format (0-based)
bed = coding_genes[['seqname', 'start', 'end', 'gene_name', 'score', 'strand']].copy()
bed['start'] = bed['start'] - 1
bed['score'] = 0
bed.to_csv('coding_genes.bed', sep='\t', header=False, index=False)
```

## Common Tasks

### Get Promoters

```python
import gtfparse
import pandas as pd

df = gtfparse.read_gtf('annotation.gtf')
tss = df[df['feature'] == 'transcript'][['seqname', 'start', 'end', 'strand', 'gene_name']]

# Get TSS position
tss['tss'] = tss.apply(lambda x: x['start'] if x['strand'] == '+' else x['end'], axis=1)

# Create promoter regions (2kb upstream)
tss['prom_start'] = tss.apply(
    lambda x: x['tss'] - 2001 if x['strand'] == '+' else x['tss'] - 1, axis=1
)
tss['prom_end'] = tss.apply(
    lambda x: x['tss'] - 1 if x['strand'] == '+' else x['tss'] + 2000, axis=1
)

promoters = tss[['seqname', 'prom_start', 'prom_end', 'gene_name']].drop_duplicates()
promoters.to_csv('promoters.bed', sep='\t', header=False, index=False)
```

### Get Introns

```bash
# Using bedtools
# First get gene body, then subtract exons
awk '$3 == "gene"' annotation.gtf | cut -f1,4,5,7 > genes.bed
awk '$3 == "exon"' annotation.gtf | cut -f1,4,5,7 > exons.bed
bedtools subtract -a genes.bed -b exons.bed > introns.bed
```

## Common Issues

### Memory Error with Large GTF

**Problem:** gtfparse loads entire file into memory
**Solutions:**
1. Filter features during load: `read_gtf('file.gtf', features=['gene'])`
2. Use gffutils database for large files

### Chromosome Naming

**Problem:** chr1 vs 1 mismatch
**Solution:** Standardize names

```python
# Add 'chr' prefix
df['seqname'] = 'chr' + df['seqname'].astype(str)

# Or remove 'chr' prefix
df['seqname'] = df['seqname'].str.replace('chr', '')
```

### GTF Parsing Errors

**Problem:** Malformed GTF attributes
**Solution:** Validate first

```bash
gffread -E annotation.gtf 2>&1 | head
```

## Resources

- [GENCODE](https://www.gencodegenes.org/) - Human/mouse annotations
- [Ensembl FTP](https://ftp.ensembl.org/) - Multi-species annotations
- [gffread documentation](http://ccb.jhu.edu/software/stringtie/gff.shtml)
- [gtfparse GitHub](https://github.com/openvax/gtfparse)
- [gffutils documentation](https://daler.github.io/gffutils/)

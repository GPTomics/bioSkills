# GTF/GFF Handling - Usage Guide

## Overview
GTF (Gene Transfer Format) and GFF3 (General Feature Format) are standard formats for gene annotations. This skill covers parsing, querying, and converting these files using gtfparse, gffutils, and gffread.

## Prerequisites
```bash
# gffread (CLI - format conversion)
conda install -c bioconda gffread

# Python packages
pip install gtfparse gffutils
```

## Quick Start
Tell your AI agent what you want to do:
- "Extract all protein-coding gene coordinates from my GTF file"
- "Get promoter regions for all transcripts"
- "Convert my GTF file to GFF3 format"

## Example Prompts

### Extracting Features
> "Extract all protein-coding genes from gencode.v44.annotation.gtf"
> "Get all exons for gene TP53 from my GTF file"
> "List all transcript IDs for a specific gene"

### Creating Derived Regions
> "Create a BED file of promoters (2kb upstream of TSS) from my GTF"
> "Extract intron coordinates from my annotation file"
> "Get 5' UTR regions for all transcripts"

### Format Conversion
> "Convert my GENCODE GTF to GFF3 format"
> "Convert gene annotations from GTF to BED format"

### Filtering and Querying
> "Filter my GTF to keep only genes on chromosome 1"
> "Get all genes with biotype 'lncRNA'"
> "Find all overlapping genes in my annotation"

## What the Agent Will Do
1. Parse the GTF/GFF file using appropriate library
2. Filter or query based on your criteria
3. Transform coordinates if needed (e.g., 1-based to 0-based for BED)
4. Output in requested format

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

## Tips
- Use `gtfparse.read_gtf('file.gtf', features=['gene'])` to filter during load and save memory
- For very large files, use gffutils to create a queryable database
- Always check chromosome naming when combining with other data (chr1 vs 1)
- Validate GTF files with `gffread -E annotation.gtf` before processing
- Remember to subtract 1 from start when converting to BED format

## Resources
- [GENCODE](https://www.gencodegenes.org/) - Human/mouse annotations
- [Ensembl FTP](https://ftp.ensembl.org/) - Multi-species annotations
- [gffread documentation](http://ccb.jhu.edu/software/stringtie/gff.shtml)
- [gtfparse GitHub](https://github.com/openvax/gtfparse)
- [gffutils documentation](https://daler.github.io/gffutils/)

# UMI Processing - Usage Guide

## Overview
Unique Molecular Identifiers (UMIs) are random sequences added to molecules before PCR amplification. They enable distinguishing PCR duplicates from biological duplicates, crucial for accurate quantification in RNA-seq, targeted sequencing, and single-cell applications.

## Prerequisites
```bash
conda install -c bioconda umi_tools
```

## Quick Start
Tell your AI agent what you want to do:
- "Extract UMIs from my reads and deduplicate after alignment"
- "Process my UMI-tagged library for accurate quantification"
- "Remove PCR duplicates using UMI information"

## Example Prompts

### Standard Workflow
> "Extract 8bp UMIs from read 1 and move them to the read header"

> "Deduplicate my aligned BAM file using UMI information"

### Library-Specific
> "Process my 10X Genomics library with cell barcodes and UMIs"

> "Handle my NEBNext library with 8bp UMIs"

### Deduplication Options
> "Use directional deduplication method for my RNA-seq data"

> "Deduplicate with cluster method for high-error-rate UMIs"

## What the Agent Will Do
1. Extract UMIs from reads and add to read headers
2. Pass through alignment with UMI information preserved
3. Deduplicate aligned reads based on UMI + alignment position
4. Generate deduplication statistics
5. Output deduplicated BAM for downstream analysis

## UMI Pattern Syntax

| Pattern | Description |
|---------|-------------|
| `N` | UMI base |
| `C` | Cell barcode base |
| `X` | Base to discard |
| `NNNNNNNN` | 8bp UMI |
| `CCCCCCCCNNNNNNNN` | Cell barcode + UMI |

## Common Library Types

| Library | Pattern |
|---------|---------|
| NEBNext | `NNNNNNNN` (8bp in R1) |
| 10X 3' v3 | `CCCCCCCCCCCCCCCCNNNNNNNNNNNN` (16bp CB + 12bp UMI) |
| TruSeq UMI | `NNNNNNNNN` (9bp in index) |

## Tips
- Extract UMIs before alignment; they must be in read headers for deduplication
- Directional deduplication is recommended for most RNA-seq applications
- High deduplication rates (>70%) may indicate library over-amplification
- Low deduplication rates (<10%) may indicate under-sequencing
- Check UMI diversity in the deduplication stats to assess library complexity

## Resources
- [umi_tools Documentation](https://umi-tools.readthedocs.io/)
- [umi_tools Publication](https://doi.org/10.1101/gr.209601.116)

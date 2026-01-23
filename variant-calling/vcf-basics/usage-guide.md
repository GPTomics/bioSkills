# VCF/BCF Basics - Usage Guide

## Overview

View, query, and understand VCF/BCF variant files using bcftools and cyvcf2. Essential for working with variant calling results.

## Prerequisites

```bash
conda install -c bioconda bcftools htslib
pip install cyvcf2
```

## Quick Start

Tell your AI agent what you want to do:

- "View the first 20 variants in my VCF file"
- "List all samples in this VCF"
- "Extract chromosome, position, and genotypes to a TSV"
- "Convert my VCF to BCF format"
- "Query variants in the region chr1:1000000-2000000"

## Example Prompts

### Viewing Files

> "Show me the header of my VCF file"

> "View variants in chromosome 1 between positions 1M and 2M"

> "How many SNPs and indels are in this VCF?"

### Extracting Data

> "Extract variant positions and allele frequencies to a table"

> "Get genotypes for all samples as a TSV file"

> "List all variant IDs (rsIDs) in the file"

### Format Conversion

> "Compress my VCF with bgzip and index it"

> "Convert this VCF to BCF format for faster processing"

### Python Analysis

> "Iterate through variants and count heterozygous calls per sample"

> "Filter variants with QUAL > 30 and write to a new file"

## What the Agent Will Do

1. Identify the VCF/BCF file location and check if it's indexed
2. Select appropriate tool (bcftools for CLI, cyvcf2 for Python)
3. Construct the query or view command with proper flags
4. Execute the operation and format output as requested
5. For region queries, ensure the file is indexed first

## Tips

- Always use bgzip (not gzip) for VCF compression - bcftools requires it
- Index files (.tbi or .csi) enable fast region queries
- BCF format is faster to process than VCF for large files
- Use `-H` flag with bcftools query to skip the header line
- cyvcf2 is faster than PyVCF for large files
- Filter status of `None` in cyvcf2 means the variant passed all filters

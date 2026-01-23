# ClinVar Lookup - Usage Guide

## Overview

Query ClinVar for variant pathogenicity classifications, review status, and disease associations using REST API or local VCF files.

## Prerequisites

```bash
pip install requests cyvcf2 pandas

# For local ClinVar queries
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
```

## Quick Start

Tell your AI agent:
- "Look up this variant in ClinVar"
- "Find pathogenic variants in BRCA1 from ClinVar"
- "Annotate my VCF with ClinVar classifications"
- "Check ClinVar review status for these variants"

## Example Prompts

### Single Variant Queries

> "What is the ClinVar classification for BRAF V600E?"

> "Look up rs121913527 in ClinVar and show disease associations"

> "Check if chr17:7577120:C>T is pathogenic in ClinVar"

### Gene-Level Queries

> "Find all pathogenic/likely pathogenic variants in TP53 from ClinVar"

> "How many ClinVar submissions are there for CFTR?"

> "List ClinVar variants in BRCA2 with expert panel review"

### Batch Annotation

> "Annotate my VCF file with ClinVar pathogenicity"

> "Add ClinVar classifications to these 100 variants"

> "Filter my variants to keep only those with ClinVar entries"

## What the Agent Will Do

1. Determine query method (API vs local VCF based on variant count)
2. Format variant identifiers for ClinVar query
3. Execute search or lookup query
4. Parse clinical significance and review status
5. Return structured results with disease associations

## Tips

- **Use local VCF** for batch queries (>100 variants) - faster than API
- **Check review status** - 2+ stars indicates higher confidence
- **Handle conflicting interpretations** - multiple labs may disagree
- **Update ClinVar regularly** - new submissions added monthly
- **CLNSIG can be multi-valued** - parse carefully for mixed classifications

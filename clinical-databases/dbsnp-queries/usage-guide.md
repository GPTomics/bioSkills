# dbSNP Queries - Usage Guide

## Overview

Query dbSNP for rsID lookups, coordinate mapping, and variant annotations using myvariant.info or NCBI Entrez APIs.

## Prerequisites

```bash
pip install myvariant biopython requests pandas
```

## Quick Start

Tell your AI agent:
- "Look up rsID rs121913527 in dbSNP"
- "Convert these genomic coordinates to rsIDs"
- "Get dbSNP annotations for my variant list"
- "Map rsIDs to GRCh38 coordinates"

## Example Prompts

### rsID Lookups

> "What are the details for rs121913527 in dbSNP?"

> "Look up these 100 rsIDs and get their genomic coordinates"

> "Find the gene associated with rs1800566"

### Coordinate Mapping

> "What is the rsID for chr7:140453136:A>T?"

> "Convert my VCF coordinates to rsIDs"

> "Map these GRCh38 positions to dbSNP rsIDs"

### Cross-Reference

> "Get ClinVar and gnomAD info for these rsIDs"

> "Find dbSNP class (SNV, indel) for my variants"

> "Check if these rsIDs are validated in dbSNP"

## What the Agent Will Do

1. Parse variant identifier (rsID or coordinates)
2. Query myvariant.info or NCBI API
3. Extract dbSNP annotations (coordinates, gene, class)
4. Map between rsIDs and genomic coordinates if needed
5. Return structured results with cross-references

## Tips

- **myvariant.info** is faster for batch queries than NCBI Entrez
- **rsIDs can be merged** - old rsIDs redirect to current ones
- **dbSNP includes all variants** - not just clinically relevant ones
- **Use SPDI notation** for precise allele representation
- **Check both hg19 and hg38** coordinates - many rsIDs have both

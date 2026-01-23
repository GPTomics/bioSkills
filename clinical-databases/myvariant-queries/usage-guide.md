# MyVariant.info Queries - Usage Guide

## Overview

Query the myvariant.info API to retrieve aggregated variant annotations from ClinVar, gnomAD, dbSNP, COSMIC, CADD, and other databases in a single request.

## Prerequisites

```bash
pip install myvariant pandas
```

## Quick Start

Tell your AI agent:
- "Query myvariant.info for this list of variants"
- "Get ClinVar and gnomAD annotations for my variants"
- "Search for pathogenic BRCA1 variants in myvariant"
- "Batch annotate these rsIDs with clinical data"

## Example Prompts

### Single Variant Queries

> "Get all annotations for rs121913527 from myvariant.info"

> "Query the BRAF V600E variant and show ClinVar classification"

> "Look up chr7:g.140453136A>T in myvariant"

### Batch Queries

> "Annotate these 500 variants with ClinVar, gnomAD, and CADD scores"

> "Query myvariant for all variants in my VCF and extract pathogenicity"

> "Batch query these rsIDs and create a summary table"

### Search Queries

> "Find all pathogenic variants in TP53 using myvariant"

> "Search for variants in the BRCA1 gene region"

> "Find COSMIC mutations in KRAS"

## What the Agent Will Do

1. Initialize myvariant client
2. Format variant identifiers (HGVS, rsID, or gene:protein)
3. Execute single or batch queries with appropriate fields
4. Parse nested JSON responses to extract key annotations
5. Handle missing data and format results as DataFrame

## Tips

- **Use HGVS notation** for most reliable matching (chr7:g.140453136A>T)
- **Batch queries** support up to 1000 variants per request
- **Specify fields** to reduce response size and improve speed
- **Check for None** values when accessing nested fields
- **myvariant.info** aggregates 40+ data sources - specify which you need

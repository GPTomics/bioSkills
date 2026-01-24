---
name: bio-clinical-databases-myvariant-queries
description: Query myvariant.info API for aggregated variant annotations from multiple databases (ClinVar, gnomAD, dbSNP, COSMIC, etc.) in a single request. Use when annotating variants with clinical and population data from multiple sources simultaneously.
tool_type: python
primary_tool: myvariant
---

# MyVariant.info Queries

## Required Imports

```python
import myvariant
```

## Initialize Client

```python
mv = myvariant.MyVariantInfo()
```

## Query Single Variant

```python
# Query by HGVS notation (recommended)
result = mv.getvariant('chr7:g.140453136A>T')

# Query by rsID
result = mv.getvariant('rs121913527')

# Query by gene and protein change
result = mv.getvariant('BRAF:p.V600E')
```

## Query Multiple Variants

```python
variants = [
    'chr7:g.140453136A>T',
    'chr17:g.7577120C>T',
    'rs121913527'
]

# Batch query (up to 1000 variants per request)
results = mv.getvariants(variants)

# With specific fields
results = mv.getvariants(
    variants,
    fields=['clinvar', 'gnomad_exome', 'dbsnp']
)
```

## Search Variants

```python
# Search by gene
results = mv.query('clinvar.gene.symbol:BRCA1', size=100)

# Search pathogenic variants in gene
results = mv.query(
    'clinvar.gene.symbol:BRCA1 AND clinvar.clinical_significance:Pathogenic',
    size=100
)

# Search by genomic region
results = mv.query('chr7:140400000-140500000')
```

## Available Fields

Common field paths for annotations:

| Field | Description |
|-------|-------------|
| `clinvar` | ClinVar annotations |
| `gnomad_exome` | gnomAD exome frequencies |
| `gnomad_genome` | gnomAD genome frequencies |
| `dbsnp` | dbSNP annotations |
| `cosmic` | COSMIC cancer mutations |
| `cadd` | CADD deleteriousness scores |
| `dbnsfp` | dbNSFP functional predictions |
| `snpeff` | SnpEff annotations |

## Extract Specific Annotations

```python
result = mv.getvariant('chr7:g.140453136A>T')

# ClinVar classification
clinvar_sig = result.get('clinvar', {}).get('clinical_significance')

# gnomAD allele frequency
gnomad_af = result.get('gnomad_exome', {}).get('af', {}).get('af')

# CADD score
cadd_phred = result.get('cadd', {}).get('phred')
```

## Batch Processing with DataFrame

```python
import pandas as pd

variants = ['rs121913527', 'rs1800566', 'rs104894155']
results = mv.getvariants(variants, fields=['clinvar', 'gnomad_exome'])

records = []
for r in results:
    records.append({
        'query': r.get('query'),
        'clinvar_sig': r.get('clinvar', {}).get('clinical_significance'),
        'gnomad_af': r.get('gnomad_exome', {}).get('af', {}).get('af')
    })

df = pd.DataFrame(records)
```

## Rate Limiting

```python
# myvariant handles rate limiting automatically
# For large batches, use chunks
def batch_query(variants, chunk_size=1000):
    all_results = []
    for i in range(0, len(variants), chunk_size):
        chunk = variants[i:i + chunk_size]
        results = mv.getvariants(chunk)
        all_results.extend(results)
    return all_results
```

## Related Skills

- clinvar-lookup - Detailed ClinVar queries
- gnomad-frequencies - gnomAD-specific frequency queries
- dbsnp-queries - dbSNP rsID lookups

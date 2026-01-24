---
name: bio-clinical-databases-clinvar-lookup
description: Query ClinVar for variant pathogenicity classifications, review status, and disease associations via REST API or local VCF. Use when determining clinical significance of variants for diagnostic or research purposes.
tool_type: python
primary_tool: requests
---

# ClinVar Lookup

## REST API Queries

### Query by Variant ID

```python
import requests

def query_clinvar_by_id(variation_id):
    '''Query ClinVar by variation ID'''
    url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi'
    params = {
        'db': 'clinvar',
        'id': variation_id,
        'retmode': 'json'
    }
    response = requests.get(url, params=params)
    return response.json()

result = query_clinvar_by_id('16609')
```

### Search by Gene

```python
def search_clinvar_gene(gene_symbol, pathogenic_only=False):
    '''Search ClinVar for variants in a gene'''
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'

    term = f'{gene_symbol}[gene]'
    if pathogenic_only:
        term += ' AND pathogenic[clinical_significance]'

    params = {
        'db': 'clinvar',
        'term': term,
        'retmax': 500,
        'retmode': 'json'
    }
    response = requests.get(url, params=params)
    return response.json()
```

### Search by HGVS

```python
def search_clinvar_hgvs(hgvs):
    '''Search ClinVar by HGVS notation'''
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
    params = {
        'db': 'clinvar',
        'term': f'{hgvs}[variant name]',
        'retmode': 'json'
    }
    response = requests.get(url, params=params)
    return response.json()
```

## Local ClinVar VCF

### Download ClinVar VCF

```bash
# GRCh38
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi

# GRCh37
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz
```

### Query Local ClinVar with cyvcf2

```python
from cyvcf2 import VCF

clinvar = VCF('clinvar.vcf.gz')

def lookup_variant(chrom, pos, ref, alt):
    '''Look up variant in local ClinVar VCF'''
    region = f'{chrom}:{pos}-{pos}'
    for variant in clinvar(region):
        if variant.REF == ref and alt in variant.ALT:
            return {
                'clnsig': variant.INFO.get('CLNSIG'),
                'clnrevstat': variant.INFO.get('CLNREVSTAT'),
                'clndn': variant.INFO.get('CLNDN'),
                'clnvc': variant.INFO.get('CLNVC')
            }
    return None

result = lookup_variant('7', 140453136, 'A', 'T')
```

## Clinical Significance Categories

| Value | Interpretation |
|-------|----------------|
| Pathogenic | Disease-causing |
| Likely_pathogenic | Probably disease-causing |
| Uncertain_significance | VUS - unknown |
| Likely_benign | Probably not disease-causing |
| Benign | Not disease-causing |
| Conflicting_interpretations | Multiple labs disagree |

## Review Status Stars

| Stars | Review Status |
|-------|---------------|
| 4 | Practice guideline |
| 3 | Expert panel reviewed |
| 2 | Multiple submitters, criteria provided |
| 1 | Single submitter, criteria provided |
| 0 | No assertion criteria |

## Parse ClinVar INFO Fields

```python
def parse_clinvar_significance(clnsig):
    '''Parse ClinVar CLNSIG field'''
    pathogenic_terms = ['Pathogenic', 'Likely_pathogenic']
    benign_terms = ['Benign', 'Likely_benign']

    if any(term in clnsig for term in pathogenic_terms):
        return 'pathogenic'
    elif any(term in clnsig for term in benign_terms):
        return 'benign'
    elif 'Conflicting' in clnsig:
        return 'conflicting'
    else:
        return 'vus'
```

## Batch Annotation with bcftools

```bash
# Annotate VCF with ClinVar
bcftools annotate \
    -a clinvar.vcf.gz \
    -c INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNDN \
    input.vcf.gz \
    -o annotated.vcf.gz
```

## Related Skills

- myvariant-queries - Aggregated queries including ClinVar
- variant-prioritization - Filter by ClinVar significance
- variant-calling/clinical-interpretation - ACMG guidelines

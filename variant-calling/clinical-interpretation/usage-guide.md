# Clinical Interpretation Usage Guide

## Overview

Prioritize and classify variants for clinical significance using ClinVar annotations, population frequencies, and ACMG classification guidelines.

## Prerequisites

```bash
# InterVar for ACMG classification
pip install intervar

# ClinVar database
wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
```

## Interpretation Framework

```
Annotated VCF
    ├── ClinVar lookup (known clinical assertions)
    ├── Population frequency (gnomAD)
    ├── Computational predictions (SIFT, PolyPhen, CADD)
    ├── ACMG classification
    └── Prioritized variant list
```

## Quick Start

### Add ClinVar Annotations

```bash
bcftools annotate \
    -a clinvar.vcf.gz \
    -c INFO/CLNSIG,INFO/CLNDN \
    input.vcf.gz -Oz -o with_clinvar.vcf.gz
```

### Filter Clinically Relevant

```python
from cyvcf2 import VCF

for v in VCF('annotated.vcf.gz'):
    clnsig = v.INFO.get('CLNSIG', '')
    if 'Pathogenic' in clnsig or 'Likely_pathogenic' in clnsig:
        print(f'{v.CHROM}:{v.POS} {v.REF}>{v.ALT[0]} - {clnsig}')
```

## ACMG Classification

The 5-tier system for variant pathogenicity:

| Class | Meaning | Action |
|-------|---------|--------|
| Pathogenic | Disease-causing | Report, clinical action |
| Likely pathogenic | Probably disease-causing | Report, consider action |
| VUS | Uncertain significance | Report with caution |
| Likely benign | Probably not disease-causing | May report |
| Benign | Not disease-causing | Usually not reported |

## Prioritization Criteria

```python
def prioritize_variant(v):
    score = 0

    # ClinVar pathogenic
    if 'Pathogenic' in v.INFO.get('CLNSIG', ''):
        score += 10

    # Rare in gnomAD
    af = v.INFO.get('gnomAD_AF', 1.0)
    if af < 0.001:
        score += 5

    # High CADD score
    cadd = v.INFO.get('CADD_PHRED', 0)
    if cadd > 20:
        score += 3

    # Protein-altering
    if 'missense' in v.INFO.get('Consequence', ''):
        score += 2

    return score
```

## Key Databases

| Database | Purpose |
|----------|---------|
| ClinVar | Clinical variant assertions |
| gnomAD | Population allele frequencies |
| OMIM | Disease-gene associations |
| CADD | Deleteriousness scores |

## Reporting Considerations

- Always include variant evidence and classification rationale
- VUS variants may be reclassified with new evidence
- Consider family segregation data when available
- Follow laboratory/institutional reporting guidelines

## Example Prompts

> "Annotate my VCF with ClinVar clinical significance"

> "Filter for pathogenic and likely pathogenic variants"

> "Prioritize variants by clinical relevance"

> "Apply ACMG classification criteria to my variants"

## What the Agent Will Do

1. Annotate variants with ClinVar clinical significance and disease associations
2. Filter by population frequency using gnomAD to flag rare variants
3. Apply computational prediction scores (CADD, SIFT, PolyPhen) for prioritization
4. Classify variants using ACMG 5-tier criteria where applicable
5. Generate a prioritized variant list ranked by clinical relevance

## Tips

- Always include variant evidence and classification rationale
- VUS variants may be reclassified with new evidence -- periodic re-analysis is recommended
- Consider family segregation data when available
- Follow laboratory/institutional reporting guidelines
- ClinVar assertions vary in review status -- prioritize variants with multiple submitters or expert panel review

## Related Skills

- variant-calling/variant-annotation - Add annotations from ClinVar, gnomAD, and other databases
- variant-calling/filtering-best-practices - Filter variants before clinical review
- variant-calling/variant-normalization - Normalize variants before database lookups
- pathway-analysis/go-enrichment - Functional enrichment of genes harboring clinical variants

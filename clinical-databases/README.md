# clinical-databases

## Overview

Query clinical and population genetics databases including ClinVar, dbSNP, gnomAD, and OMIM for variant interpretation and prioritization.

**Tool type:** python | **Primary tools:** myvariant, requests, pandas

## Skills

| Skill | Description |
|-------|-------------|
| myvariant-queries | Query myvariant.info API for aggregated variant annotations |
| clinvar-lookup | ClinVar pathogenicity classification and review status |
| gnomad-frequencies | Population allele frequencies from gnomAD |
| dbsnp-queries | rsID lookup and variant annotation |
| variant-prioritization | Filter variants by pathogenicity and frequency criteria |
| pharmacogenomics | PharmGKB/CPIC drug-gene interactions and star allele interpretation |
| polygenic-risk | PRS calculation with PRSice-2, LDpred2 |
| tumor-mutational-burden | TMB calculation with panel normalization |
| hla-typing | HLA allele calling with OptiType, HLA-HD |
| somatic-signatures | Mutational signature extraction with SigProfiler |

## Example Prompts

- "Get ClinVar classifications for my list of variants"
- "Find gnomAD allele frequencies for variants in my VCF"
- "Look up rsIDs for these genomic positions"
- "Filter my variants to keep only pathogenic/likely pathogenic"
- "Prioritize variants with AF < 0.01 and ClinVar pathogenic"
- "Get OMIM disease associations for these genes"
- "Batch annotate my VCF with ClinVar and gnomAD"
- "What is my CYP2D6 metabolizer status for *1/*4?"
- "Calculate polygenic risk score from my GWAS summary stats"
- "Calculate tumor mutational burden from my somatic VCF"
- "Run HLA typing on my WES data"
- "Extract mutational signatures from my tumor samples"

## Requirements

```bash
pip install myvariant requests pandas cyvcf2 SigProfilerExtractor
```

```r
# For LDpred2
install.packages('bigsnpr')
```

## Related Skills

- **variant-calling** - Upstream variant calling and annotation
- **database-access** - General database query patterns with Entrez
- **variant-calling/clinical-interpretation** - ACMG classification guidelines

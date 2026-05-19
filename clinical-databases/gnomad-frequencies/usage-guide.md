# gnomAD Frequencies - Usage Guide

## Overview

Query gnomAD v4 (807k samples), v3.1.2 genomes, and v2.1.1 for allele frequencies, grpmax FAF95, LOEUF gene constraint, SV/CNV catalogs, and mtDNA frequencies. The skill emphasizes the Whiffin 2017 max-credible-AF framework, the v2/v3/v4 use-case decision matrix, bottleneck-group exclusion from grpmax, and the LOEUF v2 vs v4 distribution shift.

## Prerequisites

```bash
pip install requests myvariant pandas

# Hail for bulk queries (optional)
pip install hail
gcloud auth application-default login
```

## Quick Start

Tell the agent what to do:
- "What is the v4 grpmax FAF95 for this variant? Apply ACMG BS1/BA1"
- "Filter my VCF to variants with grpmax FAF95 < 0.0001 using gnomAD v4 Hail Table"
- "Get LOEUF for these 50 candidate disease genes; note v4 has no chrX/Y constraint"
- "Compute Whiffin max-credible-AF for hypertrophic cardiomyopathy and apply BS1 to my variant list"
- "Compare v2 vs v4 AF for this variant and reconcile the discrepancy"

## Example Prompts

### Single Variant Frequency

> "Look up chr17:43094464:G:A in gnomAD v4. Report exome AF, grpmax FAF95, grpmax ancestry, and apply BS1 against breast-cancer-specific max-credible-AF."

> "What is the gnomAD v4 grpmax FAF95 for rs334 (HbS)? Note: this is a founder allele in AFR, so check that bottleneck groups are excluded properly."

> "For BRAF V600E, give me v4 exome AF, v2 AF, and reconcile any difference."

### Constraint Metrics

> "Pull LOEUF for these 30 candidate genes. Use v4 March 2024 release for autosomes; fall back to v2.1.1 for chrX/Y."

> "Rank my candidate gene list by LOEUF decile (first decile = strongly LoF-intolerant)."

> "What is the missense Z-score and missense O/E for SCN2A? Is it among the top 1% missense-constrained?"

### ACMG BS1/BA1 Application

> "Apply Whiffin 2017 max-credible-AF framework to these candidate variants for autosomal dominant cardiomyopathy. Disease prevalence 1/500, penetrance 0.8."

> "For each variant, report whether grpmax FAF95 triggers BA1 (>5%), BS1 (>max-credible), PM2_Supporting (absent/ultra-rare), or no criterion."

### Population-Specific

> "For variant chr2:165996624:T:C, compare AFs across all v4 ancestry groups. Flag if FIN/AMI/ASJ/REMAINING (bottleneck) inflate grpmax."

> "Find variants where v4 AFR FAF95 is >10x higher than NFE FAF95."

### Bulk Filtering

> "Using the gnomAD v4 exomes Hail Table at `gs://gcp-public-data--gnomad/release/4.1/ht/exomes/`, filter my cohort VCF to retain only variants with grpmax FAF95 < 0.0001."

> "Annotate my somatic VCF with gnomAD v4 grpmax FAF95 to identify likely germline contamination."

### SV / CNV / mtDNA

> "Look up these structural variants in gnomAD-SV v4. Report breakpoint AF and clinical relevance."

> "Find rare CNVs (AF < 1%) overlapping my candidate gene panel from gnomAD-CNV v4."

> "Pull mtDNA variant frequencies from gnomAD v3.1 for this Leigh syndrome candidate set."

## What the Agent Will Do

1. Choose the appropriate dataset (v2.1.1 / v3.1.2 / v4) based on use case: constraint needs v2; non-coding rare variants v3 or v4 genomes; rare-variant filtering v4 exomes.
2. Query via GraphQL for single variants or Hail Table for bulk; use myvariant.info as aggregator if multi-database overlay needed.
3. Extract grpmax FAF95 (not raw AF) for ACMG application.
4. Apply Whiffin max-credible-AF formula when gene-specific BS1 needed; default BA1 = 5% per ClinGen SVI.
5. Cross-check bottleneck-group inclusion (AMI/ASJ/FIN/REMAINING excluded from grpmax by design).
6. For chrX/Y constraint, fall back to v2.1.1 (v4 not released).
7. Pin VEP version (v4 = VEP 105; v2 = VEP 85) for consequence prediction reproducibility.

## Tips

- Use `grpmax_faf95`, not raw AF, for ACMG BS1/BA1; this is the ClinGen SVI recommendation.
- v4 genomes are the same 76,215 v3 samples reprocessed; not independent; for true non-overlap use `non_v2` subset.
- v4 includes 416,555 UK Biobank exomes; for ancestry-balanced analysis use `non_ukb` subset.
- LOEUF first decile = strongly LoF-intolerant; threshold shifted v2 < 0.35 -> v4 < 0.6 due to larger sample. Compare deciles, not absolute values.
- v4 constraint release (March 2024) is autosomes only; chrX/Y constraint requires v2.1.1.
- gnomAD-SV v4 = GRCh38, 63,046 samples; gnomAD-CNV v4 = 464,297 exome-derived CNVs. Choose by data type.
- mtDNA frequencies exist only in v3.1.2 (Laricchia 2022); apply non-Mendelian inheritance carefully.
- Filtering allele frequency formula: `(prevalence x heterogeneity x allelic_contribution) / (penetrance x 2)`. Gene-specific BS1.
- Bottleneck groups (AMI, ASJ, FIN, REMAINING) are excluded from grpmax to avoid founder-variant false BS1/BA1.
- v4.1 (May 2024) fixed AN under-counting in v4.0; rare-variant AFs inflated 5-10% in v4.0; always use v4.1 or later.
- VEP version pinning matters: a variant's consequence prediction can flip between v2 and v4 due to transcript-set updates.

## Related Skills

- clinical-databases/clinvar-lookup - Pathogenicity (gnomAD AF used for BS1/BA1 cross-check)
- clinical-databases/acmg-classification - Whiffin FAF95 framework in ACMG context
- clinical-databases/variant-prioritization - Rare-disease pipeline using grpmax_faf95
- clinical-databases/myvariant-queries - Aggregated queries including gnomAD overlay
- population-genetics/population-structure - Population stratification background

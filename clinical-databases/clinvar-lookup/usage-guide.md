# ClinVar Lookup - Usage Guide

## Overview

Query ClinVar for variant pathogenicity classifications, ClinGen Variant Curation Expert Panel (VCEP) curations, and germline-vs-somatic interpretations against the ACMG/AMP framework with ClinGen SVI specifications. The skill addresses the VCV/SCV/RCV identifier hierarchy, the 2024 XML schema overhaul, star-rating override semantics, and the ClinGen Allele Registry as the canonical cross-database join key.

## Prerequisites

```bash
pip install requests cyvcf2 pandas lxml

# Local VCF (weekly snapshot; pin to monthly archive for reproducibility)
mkdir -p clinvar/$(date +%Y%m); cd clinvar/$(date +%Y%m)
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi

# For reproducibility, use the monthly archive
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2025/clinvar_20251101.vcf.gz
```

## Quick Start

Tell the agent what to do:
- "Look up the ClinVar classification for BRAF V600E and report star rating and VCEP curation status"
- "Annotate my somatic VCF with germline + oncogenicity + somatic-clinical-impact classifications from ClinVar 2024 schema"
- "Resolve these variants to ClinGen Allele Registry CA IDs for cross-database joining with gnomAD"
- "Triangulate conflicting interpretations for this VUS: list all SCVs with submitter and star, identify the highest-confidence assertion"
- "Filter my list to retain only variants with star >=2 OR ClinGen VCEP curation"

## Example Prompts

### Single Variant Queries

> "What is the ClinVar VCV classification for chr17:43094464:G:A in BRCA1, and is there an ENIGMA-VCEP curation? Report germline classification, last evaluated date, and review status."

> "Look up rs121913529 in ClinVar; show all SCV-level assertions with submitter, classification, condition, and date last evaluated."

> "Resolve chr7:140453136:A:T to ClinGen Allele Registry CA ID and report the ClinVar VariationID linkage."

### Condition-Stratified Queries

> "For BRCA2 c.5946delT, give me the RCV-level classification for each condition rather than the VCV aggregate."

> "Find all ClinVar variants in TP53 with germline Pathogenic + somatic Tier I oncogenicity classifications under the 2024 tripartite schema."

### Conflict Resolution and Star Filtering

> "List all conflicting-interpretation variants in BRCA1, group by submitter pair, and identify whether the conflict is P/LP vs LB/B (clinically meaningful) or P vs LP (often immaterial)."

> "Filter my variant list to keep only star >=2 OR ClinGen VCEP-curated; flag stale assertions older than 36 months."

### Bulk Annotation

> "Annotate this exome VCF with ClinVar germline + oncogenicity + somatic-clinical-impact + review status + VCEP affiliation using bcftools."

> "For my 50k rare-variant cohort, batch-query ClinVar via myvariant.info with `fields=clinvar.review_status,clinvar.variant_id` and merge against gnomAD grpmax FAF95."

### Cross-Database Join

> "Resolve these 200 HGVS-g variants to CA IDs and then look each up in ClinVar, gnomAD v4 exomes, COSMIC, and MAVEdb in parallel."

## What the Agent Will Do

1. Parse the input identifier (rsID / HGVS / chromosomal coords / gene-protein) and resolve to a canonical form (CA ID or VCV).
2. Choose query mode by batch size: REST esummary for <10; myvariant.info for 10-1000; local VCF or bulk XML for >1000.
3. Pull VCV-level aggregate AND condition-stratified RCV-level classifications when the user references a specific phenotype.
4. Read the 2024 tripartite (Germline / SomaticClinicalImpact / Oncogenicity) classifications separately.
5. Apply the star-rating override hierarchy (4 > 3 > 2 > 1 > 0); flag stale assertions and conflicting interpretations.
6. Cross-check pathogenic calls against gnomAD grpmax FAF95 for BS1/BA1 reconciliation.
7. Return CA ID + VCV + RCV + per-classification labels + review status + last evaluated date.

## Tips

- Use the monthly archive (first Thursday of each month) for reproducible analyses; weekly releases are not archived.
- `CLNSIG` in `clinvar.vcf.gz` is VCV-level (variant-level aggregate); for condition-specific classification parse RCV-level XML.
- 2024 XML schema replaces `<ClinVarSet>` with `<VariationArchive>`; pipelines built before September 2024 must be re-targeted.
- Star rating 3 (VCEP) supersedes lower-star records per ClinGen FDA Recognition 2018; do not auto-aggregate by date.
- `Conflicting interpretations` is 1-star (not 2-star as sometimes reported); inspect `CLNSIGCONF` to see whether the conflict is clinically meaningful.
- Use ClinGen Allele Registry CA ID (`https://reg.clinicalgenome.org/`) for any cross-database join; ClinVar VariationID was renumbered during the 2017 schema redesign.
- ClinVar somatic classifications (`ONCDN`, `SCIDN`) were added in 2024; pre-2024 pipelines miss them silently.
- Conflict resolution is slow: only ~4% of BRCA1 missense VUS conflicts resolved despite years of effort (Harrison 2017 follow-up data).
- For pathogenicity classification logic (PVS1 decision tree, Pejaver 2022 calibrated PP3/BP4 thresholds, Tavtigian point system), defer to `clinical-databases/acmg-classification`; this skill is for querying ClinVar, not classification.

## Related Skills

- clinical-databases/acmg-classification - ACMG/AMP framework, Pejaver PP3/BP4 calibration, PVS1 decision tree
- clinical-databases/myvariant-queries - Multi-database aggregation with ClinVar overlay
- clinical-databases/variant-prioritization - Rare-disease filtering pipeline
- clinical-databases/gnomad-frequencies - Population frequency for BS1/BA1 cross-check
- variant-calling/clinical-interpretation - Clinical reporting workflow

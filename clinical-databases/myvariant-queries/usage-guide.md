# MyVariant.info Queries - Usage Guide

## Overview

Query the BioThings myvariant.info aggregator for ClinVar, gnomAD, dbSNP, dbNSFP (which contains AlphaMissense, REVEL, CADD, SpliceAI, BayesDel, and 40+ predictors), COSMIC, CIViC, and other clinical variant annotations in a single batched API call. The skill addresses dbNSFP version drift (the dominant staleness vector for in-silico predictors), Pejaver 2022 calibrated PP3/BP4 thresholds, and the reproducibility pattern of recording `_meta.src` versions alongside results.

## Prerequisites

```bash
pip install myvariant requests pandas
```

## Quick Start

Tell the agent what to do:
- "Annotate these 500 variants with ClinVar + gnomAD grpmax FAF95 + AlphaMissense + REVEL + SpliceAI in one batch"
- "Search myvariant for all ClinVar Pathogenic variants in BRCA1 with star >= 2"
- "Record `_meta.src` versions for each source alongside my results for reproducibility"
- "Find variants in chr17:43000000-44000000 with CADD phred > 25 and ClinVar P/LP"
- "Batch convert these 200 HGVS strings to canonical IDs and pull aggregated annotations"

## Example Prompts

### Batch Annotation

> "For these 100 HGVS-g variants, pull ClinVar significance, ClinVar review status, gnomAD v4 grpmax FAF95, AlphaMissense score, REVEL score, CADD phred, SpliceAI DS_max, and COSMIC IDs. Return as a pandas DataFrame and record source versions."

> "Annotate my 50k-variant cohort with ClinVar + gnomAD in batches of 1000 with rate-limit-safe sleeping; merge results into the input table."

### Lucene Search

> "Search myvariant for all variants in TP53 with ClinVar Pathogenic AND CADD phred > 25; return top 500."

> "Find variants in chr17:43044294-43125483 (BRCA1 region) with AlphaMissense score > 0.564 OR REVEL > 0.7."

> "Search for ClinVar VCEP-curated variants (review_status:reviewed_by_expert_panel) in cardiac sarcomere genes."

### Reproducibility

> "Run my annotation pipeline and record per-source versions in the output JSON sidecar: dbnsfp.version, clinvar.version, gnomad_exome.version, cadd.version."

> "Compare my last analysis (dbNSFP v4.5) to today's results; report which AlphaMissense scores changed due to dbNSFP version drift."

### Single-Source Deep Dive

> "Pull all dbNSFP fields for this variant: REVEL, BayesDel, AlphaMissense, EVE, ESM1b, PrimateAI-3D, MetaRNN, VEST4, MutPred. Apply Pejaver 2022 PP3/BP4 calibration to REVEL." (Defer Pejaver calibration logic to acmg-classification.)

> "For 100 rare missense variants, batch-query myvariant for all in-silico predictors and report the top-3 most concordant tools."

## What the Agent Will Do

1. Initialize `MyVariantInfo()` client; canonicalize input IDs (HGVS preferred over rsID for unambiguous matching).
2. For < 1000 variants: single `getvariants(list, fields=...)` call with explicit field selection.
3. For > 1000: chunk to 1000 + sleep 0.5s between chunks.
4. Always include `_meta` in fields for version tracking.
5. Use Lucene query for structured search (gene + significance + region + score thresholds).
6. Parse nested response defensively (`.get('source', {})` rather than direct key access).
7. Return DataFrame + per-source version dict; emit reproducibility sidecar.
8. Defer ACMG classification logic to `clinical-databases/acmg-classification`.

## Tips

- Use HGVS-g (e.g., `chr7:g.140453136A>T`) for unambiguous matching; rsID is multi-allelic for ~6-8% of variants.
- `_meta.src.<source>.version` is the canonical version field; record alongside results for reproducibility.
- dbNSFP version lags primary tool releases by 6-18 months; AlphaMissense entered dbNSFP 4.4 (~2024). For cutting-edge predictions query the source API directly.
- AlphaMissense developer threshold (0.564) is NOT the Pejaver 2022 PP3-calibrated threshold; do not auto-apply as PP3 evidence. ClinGen has not endorsed AlphaMissense PP3 strength as of May 2026.
- Do not stack REVEL + BayesDel + VEST4 as independent PP3 evidence; they share training data per Pejaver 2022.
- Batch limit is 1000 IDs per POST; exceeding silently truncates.
- Field path errors return None, not errors; check structure with `print(mv.getvariant(test_id))` first.
- myvariant does NOT produce ACMG calls; pair with `acmg-classification` for classification.
- For PHI-sensitive workflows, switch to OpenCRAVAT (local install); myvariant requires HTTP.
- Sample overlap exists across aggregated sources; ClinVar variants often derive from gnomAD individuals; treat as not statistically independent.

## Related Skills

- clinical-databases/clinvar-lookup - Source-level ClinVar deep queries
- clinical-databases/gnomad-frequencies - Source-level gnomAD deep queries
- clinical-databases/dbsnp-queries - Source-level rsID resolution and merge-chain
- clinical-databases/acmg-classification - ACMG with Pejaver calibrated thresholds
- clinical-databases/variant-prioritization - Pipeline using aggregated annotations

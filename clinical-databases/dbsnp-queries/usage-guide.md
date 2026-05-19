# dbSNP Queries - Usage Guide

## Overview

Resolve rsIDs, navigate Build 156 JSON-based RefSNP records, follow multi-hop merge chains through `RsMergeArch` / `SNPHistory`, and normalize variant identifiers across rsID, SPDI, HGVS-g/c/p, and VCF conventions. The skill emphasizes rsID as a *cluster identifier* rather than a unique variant identifier and treats SPDI / ClinGen Allele Registry CA ID as the canonical join key.

## Prerequisites

```bash
pip install myvariant biopython requests pandas

# Optional: bulk download for genome-wide rsID lookup (each ~5-20GB compressed)
wget ftp://ftp.ncbi.nlm.nih.gov/snp/latest_release/JSON/refsnp-chr17.json.bz2
wget ftp://ftp.ncbi.nlm.nih.gov/snp/latest_release/JSON/refsnp-merged.json.bz2
wget ftp://ftp.ncbi.nlm.nih.gov/snp/latest_release/JSON/refsnp-withdrawn.json.bz2
```

## Quick Start

Tell the agent what to do:
- "Resolve rs429358 to its current canonical rsID through any merge chain"
- "Convert chr19:44908684:T:C to SPDI canonical form and find any associated rsID"
- "For these 500 deprecated rsIDs, trace merge histories and report final rsIDs"
- "Normalize this VCF's IDs by replacing rsIDs with SPDI"
- "Check whether rs1799963 is withdrawn or merged, and report ALFA + gnomAD frequencies"

## Example Prompts

### Single rsID Resolution

> "Look up rs121913529 in Build 156 and return: current canonical rsID, GRCh38 coordinates, all alleles, gene context, and merge history."

> "Resolve rs334 (HbS sickle cell variant) and report SPDI canonical form plus ALFA frequencies by ancestry."

> "Is rs1801131 (MTHFR A1298C) multi-allelic in dbSNP? Show the full alleles list."

### Merge-Chain Resolution

> "rs1799963 was reported in a 2018 paper. Trace its merge chain to the current canonical rsID."

> "For these 100 rsIDs from a 2010 candidate-gene study, identify which are still current, which have merged, and which are withdrawn."

### Coordinate Conversion

> "Convert chr17:43091031:G:A (VCF format, GRCh38) to canonical SPDI and lookup the associated rsID."

> "Given HGVS notation NM_000059.3:c.5946delT, find the genomic SPDI on GRCh38 and the corresponding rsID."

> "For my VCF of 50k variants, replace IDs with canonical SPDI; flag multi-allelic sites where rsID is ambiguous."

### Cross-Database Frequency

> "For rs6025 (factor V Leiden), give me ALFA per-ancestry frequencies AND gnomAD v4 grpmax FAF95 and reconcile any discrepancy."

> "Identify variants where ALFA AF and gnomAD AF differ by >2x and flag for review."

### Bulk Operations

> "Stream-parse `refsnp-chr17.json.bz2` to extract all rsIDs in BRCA1 + BRCA2 regions with their ALFA frequencies."

> "From `refsnp-merged.json.bz2`, build a deprecated-rsID lookup table for our cohort's old genotyping array."

## What the Agent Will Do

1. Parse the input identifier (rsID, HGVS, VCF coords, SPDI) and detect its representation.
2. Choose query mode: Variation Services REST (full JSON, single variant); myvariant.info (batch with annotation overlay); local JSON bulk (genome-wide); E-utilities (legacy, thin summary; avoid).
3. If rsID is deprecated, follow the multi-hop merge chain to the current canonical rsID.
4. If withdrawn, return `is_withdrawn=true` with withdrawal reason.
5. For coord-to-rsID conversions, use SPDI canonical form (`/spdi/{spdi}/rsids`); for rsID-to-coord, parse `primary_snapshot_data.placements_with_allele`.
6. Flag multi-allelic clusters and require allele-match filtering before downstream use.
7. Return SPDI + canonical rsID + alleles + gene + frequency (ALFA and/or gnomAD).

## Tips

- Treat rsID as a human label, not a join key. Use SPDI or CA ID for any cross-database join.
- ~6-8% of rsIDs are multi-allelic; naive `rsID -> variant` mappings silently fail.
- Variation Services REST returns the full Build 156 JSON; E-utilities (Entrez `db=snp`) returns a thin legacy summary.
- Multi-hop merges (rs3 -> rs2 -> rs1 -> rs0) are common for old rsIDs; one-hop lookup against `RsMergeArch` is insufficient.
- Withdrawn rsIDs are in `SNPHistory.bcp.gz` (or `refsnp-withdrawn.json.bz2`), not `RsMergeArch`.
- SPDI position is 0-based half-open; HGVS is 1-based fully-closed; intentional off-by-one when converting.
- ALFA (dbSNP-embedded) covers ~1M dbGaP samples and includes array-only variants; gnomAD is deeper at rare variants but sequencing-only.
- Strand-ambiguous (A/T, C/G) variants: rsID is locus-level; check ref/alt to detect strand-flips.
- For genome-wide rsID lookups, download chromosome-partitioned JSON files (`refsnp-chr{N}.json.bz2`) rather than streaming all of dbSNP.

## Related Skills

- clinical-databases/myvariant-queries - Aggregated rsID + annotation queries
- clinical-databases/clinvar-lookup - ClinVar VariationID vs rsID linkage
- clinical-databases/gnomad-frequencies - Frequency lookups by canonical SPDI
- clinical-databases/variant-prioritization - Pipeline using normalized variant IDs
- database-access/entrez-search - General Entrez query patterns

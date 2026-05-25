# Ensembl REST Usage Guide

## Overview

Query Ensembl's REST API for gene/transcript/protein lookup, sequence retrieval, Compara orthologs, VEP variant annotation, regulatory features, and LD. Encodes: the 15 req/sec rate limit (55K/hour); version-pinned archive endpoints (`https://e110.rest.ensembl.org`) for reproducibility; symbol-vs-Ensembl-ID stability (resolve symbol to ID once, then query by ID); Ensembl divisions (vertebrates main host vs Ensembl Genomes for plants/fungi/metazoa/bacteria); and the defection rules (BioMart for >5K queries; local VEP for bulk variant annotation).

## Prerequisites

```bash
pip install requests
```

No API key required. Respect `Retry-After` on 429.

## Quick Start

- "Resolve gene symbol BRCA1 to its Ensembl Gene ID; use the ID for all downstream queries"
- "Get protein sequence for ENSG00000139618 (BRCA2)"
- "What genes are in chr17:43000000-43200000? Use the overlap endpoint"
- "VEP annotate a single missense variant by region/allele notation"
- "Pin queries to release 110 archive (https://e110.rest.ensembl.org) for reproducibility"

## Example Prompts

### Resolve to stable ID first

> "Don't use gene symbols in downstream queries -- HGNC renames break workflows (MARCH1 -> MARCHF1 in 2020). Resolve via /lookup/symbol/human/BRCA1 once, persist the Ensembl Gene ID (ENSG...), and use ID-based endpoints from then on."

### Version pinning

> "I'm writing a paper. Pin all Ensembl REST calls to https://e110.rest.ensembl.org so re-running in 2030 returns the same gene models. Live https://rest.ensembl.org follows the current release and drifts."

### VEP for ad hoc variants

> "Annotate this single variant (chr17:41276135 T>G) with VEP. Use /vep/human/region/17:41276135-41276135:1/G. For bulk (>1K variants), download VEP and run locally -- REST rate limit makes bulk infeasible."

### Compara orthologs

> "Get all Compara orthologs of human BRCA1 across mouse and zebrafish. Use /homology/symbol/human/BRCA1?type=orthologues. Filter by target_species and ortholog_one2one for high-confidence calls."

### Overlap query

> "What genes are in chr17:43000000-43200000? Use /overlap/region/human/{region}?feature=gene. For multiple feature types, repeat with feature=transcript, feature=regulatory, etc."

### Non-vertebrate

> "For Arabidopsis or fungi, switch to rest.ensemblgenomes.org -- the main rest.ensembl.org host is vertebrates only."

## What the Agent Will Do

1. Always resolve gene symbol to Ensembl Gene ID at pipeline start; use IDs downstream.
2. Pin to an archive endpoint (`e110`, `e111`, etc.) for any reproducible analysis.
3. Respect 15 req/sec rate limit; sleep 0.07s between calls; handle 429 with Retry-After.
4. For >5,000 queries, recommend BioMart bulk export (separate skill).
5. For >1,000 variant annotations, recommend local VEP install.
6. For non-vertebrates, switch to `rest.ensemblgenomes.org` host.
7. Use ID-based endpoints over symbol-based for any persistent code.
8. Document the pinned release in pipeline metadata.

## Tips

- Symbol-based endpoints are convenient for interactive work; ID-based are for reproducible pipelines.
- Ensembl quarterly releases change gene model versions, exon coordinates, and sometimes Gene IDs themselves -- always pin a release for published work.
- `Accept: application/json` header is the canonical request; the REST API also supports XML and other formats but JSON is the default.
- For >5,000 lookups, the BioMart bulk export beats REST loops by 100x; see `biomart-queries`.
- VEP via REST is for ad hoc work; for production variant annotation, install VEP locally.
- The species name accepted is the Ensembl species name (`homo_sapiens` or `human`); see `/info/species` to enumerate.
- For HGNC-renamed genes (MARCH1, SEPT*), Ensembl mirrors the rename. Symbol lookup for the old name returns 404; use the new symbol or pre-resolve.
- Some endpoints accept multiple IDs in POST body (e.g. `/lookup/id` POST with `{"ids": [...]}`); use POST for batches of <500 IDs to reduce request count.
- For older assembly versions (GRCh37), use `https://grch37.rest.ensembl.org`.

## Related Skills

- biomart-queries - Bulk ID export via BioMart (>5K queries)
- ortholog-inference - Compara orthologs and other ortholog resources
- uniprot-access - Cross-reference Ensembl IDs in UniProt entries
- variant-calling/variant-annotation - Local VEP for bulk variant work
- ncbi-datasets-cli - NCBI alternative for genome / gene metadata
- entrez-search - NCBI alternative for non-Ensembl-native queries

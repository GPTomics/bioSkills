# Entrez Search Usage Guide

## Overview

Search NCBI databases with Biopython's `Bio.Entrez` (ESearch, EInfo, EGQuery, ESpell). The skill emphasizes the parts of the API that bite real users: the Entrez Query Translator rewriting unqualified queries, the 9,999 silent cap on non-history responses, history-server semantics (8-hour TTL, idle eviction, QueryKey chaining), `[Organism]` taxonomy walks, and the weekly index lag that breaks 'I just submitted my data' expectations.

## Prerequisites

```bash
pip install biopython
```

Set the required Entrez identifiers globally:

```python
from Bio import Entrez
Entrez.email = 'researcher@institution.edu'
Entrez.api_key = 'optional_api_key'  # raises rate from 3 to 10 req/sec
```

Get an API key at `https://www.ncbi.nlm.nih.gov/account/settings/`. Without it, the 3 req/sec ceiling makes any non-trivial workflow painful.

## Quick Start

- "Find PubMed articles about CRISPR published in 2024, return the first 50 PMIDs"
- "How many SRA runs exist for BioProject PRJNA123456?"
- "Show me which NCBI databases contain records mentioning BRCA1"
- "Search nucleotide for human RefSeq mRNAs between 500 and 5000 nt"
- "List the searchable fields for the ClinVar database"

## Example Prompts

### Building a reproducible query

> "Search PubMed for 2024 papers on tumor-mutational-burden in non-small-cell lung cancer. Print the QueryTranslation so I can lock the exact field-qualified rewrite into my code, then return the count and first 100 PMIDs."

### Avoiding the 9999 cap

> "Find every RefSeq mRNA for Homo sapiens. The total is large -- use the history server (usehistory='y') so downstream EFetch can pull in chunks without re-sending IDs."

### Cross-database discovery

> "I have the gene symbol DDX3X. Use EGQuery to show which NCBI databases contain records mentioning it, then drill into the gene database with a field-qualified search to get the canonical Gene UID."

### Chaining queries on the history server

> "Find PubMed records about 'BRCA1 AND breast cancer', then on the same WebEnv search 'review[PT] AND 2024[PDAT]'. Intersect QueryKey #1 AND #2 to get 2024 reviews about BRCA1 in breast cancer."

### Diagnosing a 'wrong count' bug

> "My esearch for 'MARCH1 AND human' returned no hits. Print the QueryTranslation, then re-run with the proper field-qualified form using the HGNC-permanent symbol MARCHF1."

## What the Agent Will Do

1. Set `Entrez.email` and (if available) `Entrez.api_key` before any call.
2. Translate natural-language requests into field-qualified Entrez query syntax.
3. Run ESearch with `retmax=0` first to get the count, then decide between direct retrieval and the history server based on the count.
4. Print `QueryTranslation` so the user can verify the EQT did the right thing.
5. For ambiguous gene symbols, look up the canonical HGNC name first.
6. Respect rate limits (0.34s no key, 0.10s with key) and surface warnings if hitting the 9999 cap.
7. Hand off WebEnv/QueryKey to downstream EFetch or ESummary workflows.

## Tips

- For multi-token or ambiguous gene symbols (MARCHF1, SEPTIN9, the Excel-renamed set), always go through `gene` db lookup first and use the resulting GeneID downstream -- gene symbols are unstable.
- `[Organism]` defaults to a taxonomy walk; for a single species use a binomial like `Homo sapiens[ORGN]`; to disable the walk use `[Organism:exp]`.
- Fresh deposits (< 48h) often aren't searchable yet but are EFetch-able by accession -- use EFetch directly when the accession is known.
- For very large result sets, push to the history server once with ESearch then iterate EFetch -- never re-send 50,000 IDs in successive calls.
- PMC subset filtering uses `pubmed pmc[sb]` not a separate db; the underlying records are PubMed UIDs.
- EGQuery counts can lag the per-database indexes by 1-2 days; for authoritative counts loop ESearch over the db list.

## Related Skills

- entrez-fetch - Retrieve actual records after searching
- entrez-link - Cross-database references (gene to protein, sequence to PubMed)
- batch-downloads - Pull large history-server result sets efficiently
- ncbi-datasets-cli - Modern CLI for genome and gene metadata
- geo-data - Specialized search for the gds (GEO) database
- sra-data - Search and download SRA runs

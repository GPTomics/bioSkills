# UniProt Access Usage Guide

## Overview

Query UniProt's REST API (`https://rest.uniprot.org/`). Encodes the 2022 endpoint migration (legacy `https://www.uniprot.org/uniprot/...` is deprecated; JSON schema changed substantially), the four UniProt resources (UniProtKB / UniRef / UniParc / Proteomes), search vs stream vs accessions endpoint trade-offs, isoform handling (`P04637-2` suffix), async ID-mapping job pattern, the critical `?fields=` for bulk pulls, and the Swiss-Prot vs TrEMBL distinction.

## Prerequisites

```bash
pip install requests pandas
```

No API key required. Rate limit generous (~200 req/sec). ID-mapping has its own job queue.

## Quick Start

- "Fetch UniProt P04637 (TP53) as JSON; extract gene name, sequence, PDB cross-refs"
- "Search for human reviewed kinases with PDB structures; return as TSV with selected fields"
- "Stream all human Swiss-Prot entries (~20K) with the `fields=` parameter to keep payload small"
- "Map 100 Ensembl Gene IDs to UniProt accessions via the async ID-mapping job"
- "Download the human reference proteome FASTA (UP000005640)"

## Example Prompts

### Single entry with defensive parsing

> "Fetch UniProt P04637 as JSON and pull: primaryAccession, uniProtkbId, proteinDescription.recommendedName.fullName.value, genes[0].geneName.value, sequence.value, and all PDB cross-refs. Use .get() chains for every nested field -- the schema is deeply nested and some fields are optional."

### Bulk pull with explicit fields

> "Search human reviewed kinases. Use /search with format=tsv and fields=accession,gene_primary,protein_name,length,xref_pdb,xref_alphafolddb,size=500. Don't omit fields= -- default returns full JSON entries (~20-30 KB each) and 500 of them is 15 MB."

### Stream for >500 results

> "I want all human Swiss-Prot entries (~20,000). The /search endpoint caps at 500/page. Use /stream with fields=accession,gene_primary,length to stream all results in one HTTP response."

### ID mapping async pattern

> "Convert 100 Ensembl Gene IDs to UniProt accessions. Use POST /idmapping/run with from=Ensembl, to=UniProtKB. Poll /idmapping/status/{jobId} every 3 seconds with a 600-second timeout (jobs can hang). Retrieve via /idmapping/results/{jobId}."

### Isoform handling

> "TP53 (P04637) has multiple isoforms. The bare P04637 fetch returns the canonical sequence only. Walk comments[type=ALTERNATIVE PRODUCTS].isoforms; for each isoform, fetch P04637-2.fasta, P04637-3.fasta, etc."

### Resolve obsolete accessions

> "These accessions from a 2015 paper may be obsolete or merged. Run ID mapping from=UniProtKB_AC-ID, to=UniProtKB to resolve to current primary accessions before fetching."

## What the Agent Will Do

1. Always use the post-2022 endpoint (`rest.uniprot.org`); flag legacy URLs as broken.
2. Pick the resource: UniProtKB (curated/auto entries), UniRef (clusters), UniParc (archival), Proteomes (organism sets).
3. For any bulk pull, specify `fields=` explicitly to control payload.
4. For >500 results, use `/stream` not `/search`.
5. For ID mapping, follow the async pattern: submit, poll with timeout, retrieve.
6. Filter `reviewed:true` for reference-quality analyses (avoids 250M TrEMBL).
7. Combine queries with `organism_id:9606` to scope by species.
8. Navigate nested JSON with `.get()` chains and defensive defaults.
9. For isoforms, fetch the `-N` suffixed accession separately; canonical is the default.

## Tips

- The 2022 endpoint migration is the biggest gotcha for legacy code. Old `https://www.uniprot.org/uniprot/...` URLs may redirect but the JSON schema is the new one -- old field paths break.
- `?fields=` is essential for bulk; without it each entry is 20-30 KB JSON.
- Swiss-Prot (`reviewed:true`) is ~570K entries; TrEMBL (auto-annotated) is ~250M. Always filter for reference work.
- The /stream endpoint has no 500-result cap and is the right call for any bulk pull.
- ID mapping is async: submit, poll, retrieve. Set a poll timeout -- the API doesn't fail-soft.
- Isoforms have `-N` suffix; the bare accession is the canonical sequence only.
- Gene symbol `gene:TP53` matches across species; combine with `organism_id:9606` for specificity.
- `gene_exact:` avoids wildcard matches that `gene:` allows.
- For reproducible analyses, pin a UniProt release in the citation (e.g. UniProt 2024_06).
- Reference proteomes (one per species) are at `/proteomes/{upid}.fasta.gz` -- e.g. UP000005640 is human.

## Related Skills

- entrez-fetch - NCBI protein records (RefSeq, GenPept) alternative
- biomart-queries - Alternative bulk ID-mapping via BioMart
- ortholog-inference - OMA needs UniProt accession resolution
- structural-biology/structure-io - Download PDB structures from UniProt xrefs
- structural-biology/alphafold-predictions - AlphaFoldDB IDs in UniProt xrefs
- pathway-analysis/go-enrichment - GO terms pulled from UniProt

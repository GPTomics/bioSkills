# Entrez Link Usage Guide

## Overview

Navigate between NCBI databases using `Bio.Entrez.elink()`. The skill encodes the most consequential decision in ELink work -- which `linkname` to use for each (`dbfrom`, `db`) pair -- and the auxiliary decisions: which `cmd` variant (`neighbor`, `neighbor_history`, `acheck`, `neighbor_score`), how to handle batches >200 IDs via EPost + history, and how to guard against the asymmetry of link tables.

## Prerequisites

```bash
pip install biopython
```

```python
from Bio import Entrez
Entrez.email = 'researcher@institution.edu'
Entrez.api_key = 'optional_api_key'  # 3 -> 10 req/sec
```

## Quick Start

- "Find RefSeq proteins for human BRCA1 (Gene UID 672) using the curated link, not the noisy all-protein link"
- "Get GeneRIF-curated genes mentioned in PubMed paper 35412348"
- "Discover all available link tables for a gene record"
- "Link 5,000 PubMed UIDs to genes -- the batch is too big for a comma-joined URL"
- "Find SRA runs for BioProject PRJNA123456"

## Example Prompts

### Picking the right linkname

> "For gene UID 672, get the linked proteins using linkname='gene_protein_refseq' so we get the 1-5 curated isoforms instead of the 500-protein all-proteins variant."

### Enumerating link options

> "Before I build this pipeline, show me every available linkname for (dbfrom=gene, source-id=672) using cmd='acheck'. I want to see what NCBI exposes and what each link's curation level is."

### Asymmetric round-trip awareness

> "Find genes mentioned in PMID 35412348 using pubmed_gene_rif (curated only). Then for each gene, find back-links to PubMed using gene_pubmed -- and warn me when the round-trip set isn't the original PMID, because pubmed_gene is text-mined and gene_pubmed is curated."

### Large batch via history server

> "Link 5,000 gene UIDs to RefSeq proteins. The id list is too long for a comma-joined URL, so EPost the UIDs in chunks of 200, then ELink with cmd='neighbor_history' and hand the resulting WebEnv to downstream EFetch."

### BioProject to SRA

> "I have BioProject PRJNA123456. Resolve to a UID, then ELink to SRA to get the run UIDs. Hand off to the sra-data skill to download FASTQ."

## What the Agent Will Do

1. Set `Entrez.email` and `Entrez.api_key`.
2. For unfamiliar (`dbfrom`, `db`) pairs, run `cmd='acheck'` first to enumerate linknames.
3. Pick a `linkname` deliberately -- prefer curated variants (`*_refseq`, `*_rif`, `*_swissprot`) for analyses; use the umbrella `gene_protein` only for exploration.
4. For >200 source IDs, EPost first and ELink with `cmd='neighbor_history'`.
5. Iterate the response as one LinkSet per input UID -- never assume a single LinkSetDb covers all inputs.
6. Guard for empty `LinkSetDb` before indexing.
7. Document the asymmetry of round-trip queries when results matter for publication.

## Tips

- The `linkname` is the most important parameter -- always specify it explicitly. The default (no linkname) returns the union of all link tables, which is rarely the desired result.
- For PubMed -> gene links, `pubmed_gene_rif` is curated and high-quality; `pubmed_gene` includes text-mining.
- For a single source ID, `record[0]['LinkSetDb'][0]['Link']` works. With multiple inputs, iterate the full `record` list -- one entry per input.
- ELink results sometimes include `Score` (for `pubmed_pubmed` `cmd='neighbor_score'`); most other link types don't.
- For very large batches (>1000 inputs), `cmd='neighbor_history'` is the only practical option -- it returns WebEnv + QueryKey for downstream EFetch.
- Asymmetric link tables are a feature, not a bug. The directional difference reflects different curation pipelines.
- `homologene` (linkname `gene_homologene`) was deprecated in 2014. Data still queryable but new entries stopped -- prefer Ensembl Compara or OrthoFinder for current orthology.

## Related Skills

- entrez-search - Resolve accessions/symbols to UIDs before linking
- entrez-fetch - Retrieve content for linked UIDs
- batch-downloads - Pull large history-server linksets efficiently
- geo-data - GEO-specific links (gds <-> sra, pubmed, bioproject)
- ncbi-datasets-cli - Modern CLI for many gene/genome cross-reference queries

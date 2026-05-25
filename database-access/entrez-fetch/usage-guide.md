# Entrez Fetch Usage Guide

## Overview

Retrieve full records (`EFetch`) or lightweight document summaries (`ESummary`) from NCBI databases using Biopython's `Bio.Entrez`. The skill encodes the rettype/retmode decision matrix per database, ESummary-vs-EFetch triage for bulk metadata work, GI deprecation (records after 2017 have only accession.version), the `gbwithparts` trap for WGS assemblies, XML schema drift, and accession-versioning for reproducible analyses.

## Prerequisites

```bash
pip install biopython
```

```python
from Bio import Entrez, SeqIO
Entrez.email = 'researcher@institution.edu'
Entrez.api_key = 'optional_api_key'  # raises rate to 10 req/sec
```

## Quick Start

- "Fetch the GenBank record for NM_007294.4 and parse the CDS coordinates"
- "Get organism + length for these 5,000 nucleotide UIDs without downloading sequence"
- "Download the CDS-translated proteins from the E. coli K-12 reference genome in one EFetch call"
- "Pull the full PubMed XML for PMID 35412348 including MeSH terms"
- "Convert a list of SRA UIDs to SRR run accessions with library size metrics"

## Example Prompts

### Choosing ESummary over EFetch

> "I have 8,000 nucleotide UIDs. I only need the organism, accession.version, and sequence length. Use ESummary in chunks of 500, not EFetch -- it's an order of magnitude cheaper."

### Reproducibility via versioned accessions

> "Fetch GenBank for NM_007294 but lock the result to whatever specific .version is current today. Save the accession.version so reruns next year fetch the same record content."

### Handling WGS records correctly

> "Fetch the genome assembly for WGS accession ABFD01000000 -- use rettype='gbwithparts' so the CONTIG sequences are inlined, otherwise the record comes back with no actual sequence."

### One-shot CDS extraction

> "Download all CDS translations from RefSeq NC_000913.3 (E. coli K-12) using rettype='fasta_cds_aa'. Don't walk the GenBank features manually -- let NCBI do the extraction server-side."

### SRA metadata conversion

> "Convert these SRA UIDs to SRR run accessions plus Bases/Spots/AvgLength metrics using EFetch with rettype='runinfo'. Parse the CSV and return as a DataFrame."

## What the Agent Will Do

1. Identify whether the request needs full record content (EFetch) or only metadata (ESummary).
2. Choose the `(rettype, retmode)` triple appropriate for the source database -- never default-guess.
3. For modern records, always pass `accession.version` instead of legacy GI numbers.
4. For WGS or assembly records, use `rettype='gbwithparts'` to inline the sequence.
5. Parse results with `SeqIO` (sequences) or `Entrez.read()` (XML); use `.get()` defensively for nested XML fields.
6. Sniff the response start (`LOCUS`, `>`, `<?xml`) before parsing -- guard against HTML error pages.
7. Respect rate limits (0.34s no key, 0.10s with key) and chunk batches to fit URL length limits.

## Tips

- ESummary is 5-10x cheaper than EFetch for metadata-only work. Always reach for it first when sequence content isn't needed.
- A bare accession `NM_007294` resolves to whatever version is current; this is fine for exploration but breaks reproducibility. Always pin `.version` for published analyses.
- For PubMed records, `rettype='medline'` parsed with `Bio.Medline.parse()` is more schema-stable than the XML route -- preferred for long-lived parsers.
- The default `gb` rettype on a WGS / assembly record returns the CONTIG join() statement but NOT the sequence. Use `gbwithparts` for assemblies; check `len(record.seq)` after parsing.
- EFetch URL has a practical ~2000 char limit; for >100 IDs in a comma-joined batch, either chunk or use the history server (push IDs via `Entrez.epost` first; see `batch-downloads`).
- For taxonomy: when given a species name, ESearch the `taxonomy` db to get the TXID, then EFetch with `db=taxonomy` for the lineage -- never assume a TXID from memory.
- XML schema drift is real. Pin BioPython in production and expect to update parsers when NCBI changes the schema (no version notice is given).

## Related Skills

- entrez-search - Find UIDs to fetch
- entrez-link - Cross-database navigation via ELink
- batch-downloads - History-server pipelines for large fetches
- ncbi-datasets-cli - Modern CLI for genome and gene metadata
- sequence-io/read-sequences - Parse downloaded FASTA/GenBank locally

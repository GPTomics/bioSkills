# Batch Downloads Usage Guide

## Overview

Bulk-download records from NCBI E-utilities efficiently. Encodes the strategy decision (direct fetch vs EPost vs history server vs Datasets CLI), precise rate-limit math, WebEnv lifecycle for long jobs, EPost 200-ID per-call limit, retry/resume design, integrity verification, and the defection rule to NCBI Datasets v2 CLI for genome/gene bulk work.

## Prerequisites

```bash
pip install biopython
# For the modern bulk path on genome/gene data:
conda install -c conda-forge ncbi-datasets-cli
```

```python
from Bio import Entrez
Entrez.email = 'researcher@institution.edu'
Entrez.api_key = 'YOUR_KEY'  # Get at ncbi.nlm.nih.gov/account/settings/; mandatory for bulk
Entrez.tool = 'project-name'
```

## Quick Start

- "Download all human RefSeq mRNAs to a single FASTA, with disk checkpointing so the job can resume"
- "Pull EFetch records for these 5,000 protein accessions; use EPost first to avoid URL-length errors"
- "Download every PubMed abstract for 'CRISPR AND 2024[PDAT]' in MEDLINE format"
- "Compare cost of pulling 100K nucleotide records via EFetch vs NCBI Datasets CLI"
- "Verify the downloaded FASTA has the expected record count and no truncation"

## Example Prompts

### Resumable bulk download

> "Download all RefSeq mRNAs for Homo sapiens to refseq_mrna.fasta. Use the history server. Checkpoint progress to ckpt.json so if my SSH session drops we resume from where we left off, not from zero."

### Large known-ID list

> "I have 5,000 protein accessions in a file. EPost them in chunks of 200 (since EPost's per-call limit is 200), then EFetch by WebEnv/QueryKey in batches of 500."

### Defecting to Datasets CLI

> "I need every RefSeq bacterial genome assembly. Don't loop EFetch -- use 'datasets download genome taxon Bacteria --refseq' from the NCBI Datasets v2 CLI. Show me the equivalent E-utils pipeline so I can see why Datasets is the right tool."

### Production retry / session expiry

> "Build a download with exponential-backoff retry on 429, detection of HTTP-200-with-ERROR-body (WebEnv expired), and automatic re-ESearch + resume from the disk checkpoint."

### Post-download integrity

> "After the download, parse the output FASTA with SeqIO and assert the record count matches what ESearch returned. If not, surface a warning."

## What the Agent Will Do

1. Pick the retrieval strategy from the decision matrix in SKILL.md based on record count, source DB, and whether IDs are known.
2. Set Entrez.email/api_key/tool before any call.
3. Use `usehistory='y'` for any ESearch expected to return >5000 records.
4. EPost ID lists >200 in chunks of 200.
5. Choose `batch_size` per rettype (500-1000 FASTA; 100-200 GB/XML).
6. Apply correct sleep (0.34s without key, 0.10s with key).
7. Implement disk checkpointing of `retstart` cursor for resumability.
8. Detect WebEnv expiry by parsing response body for `<ERROR>` (HTTP 200 is misleading).
9. On retry: jittered exponential backoff for 429; truncate-to-newline for partial-chunk recovery.
10. Defect to `datasets` CLI for genome/gene bulk work; document the choice.

## Tips

- Parallelizing API calls is the wrong bulk strategy. One stream with the history server + larger batches is faster and more polite than N parallel streams. Max ~4 workers with API key, 1 without.
- For >100,000 records of any kind, run the job outside US weekday business hours (NCBI ToS).
- WebEnv has an 8-hour absolute TTL and ~15-min idle eviction. Long jobs MUST checkpoint to disk and re-run ESearch on expiry.
- The 9,999 silent retmax cap (legacy esearch.fcgi without `usehistory='y'`) drops the rest of the result set with no error. Always use the history server above 5,000.
- For genome assemblies, NCBI Datasets v2 CLI is the supported bulk endpoint. It handles checksums and parallel download out of the box. Use E-utils only when Datasets doesn't cover the data type.
- Always verify post-download: SeqIO record count vs ESearch Count is a 1-line integrity check.
- For raw sequencing reads, batch-downloads is not the right skill -- use `sra-data` (SRA toolkit) or ENA mirror.

## Related Skills

- entrez-search - Build queries before batch-fetching
- entrez-fetch - Single-record fetches
- entrez-link - Chain ELink with neighbor_history for bulk cross-db work
- ncbi-datasets-cli - Modern bulk endpoint for genome and gene data
- sra-data - SRA toolkit for raw sequencing read downloads
- geo-data - GEO supplementary file downloads

# Read Sequences - Usage Guide

## Overview

This skill enables AI agents to help you read biological sequence data from common file formats using Biopython. It covers FASTA, FASTQ, GenBank, and 40+ other formats, and the choice between streaming, in-memory, and on-disk-indexed access.

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:

- "Read the sequences from my FASTA file"
- "Load this GenBank file and show me the gene features"
- "Count how many sequences are in this FASTQ file"
- "Get the sequence with ID 'NM_001234' from this large FASTA without loading the whole file"
- "Build a reusable index over these 10GB of FASTA files for repeated lookups"

## Example Prompts

### Basic Reading
> "Parse sequences.fasta and print each sequence ID and length"

### Working with GenBank
> "Read the GenBank file and extract all CDS product names"

### FASTQ Quality Analysis
> "Load reads.fastq and calculate the average quality score for each read"

### Large File Handling
> "I have a 10GB FASTA file and need to extract just one sequence by ID without loading the whole file"

### Reusable Index
> "I query the same set of reference FASTA files every run; build a persistent index so lookups are instant next time"

## What the Agent Will Do

1. Import Bio.SeqIO
2. Choose the access pattern: stream with `parse`, assert single with `read`, or random-access with `to_dict`/`index`/`index_db`
3. Specify the correct lowercase format string for the file type
4. Iterate or look up records as needed
5. Extract the requested information from SeqRecord objects

## Supported File Types

The agent can read these common formats:
- **FASTA** (.fasta, .fa, .fna, .faa)
- **FASTQ** (.fastq, .fq)
- **GenBank** (.gb, .gbk)
- **EMBL** (.embl)
- **Swiss-Prot** (.dat)
- **Alignment formats** (PHYLIP, Clustal, Stockholm)

## Tips

- Stream by default. `parse()` holds one record at a time and scales to any file size; only load everything into memory when random access on a small file demands it.
- Never `list()` or `to_dict()` a huge file; it can exhaust RAM. Use `index()` for random access into one large file, or `index_db()` for files larger than memory, many files at once, or an index reused across runs.
- A `parse()` generator is one-pass: it exhausts silently after the first loop and is not subscriptable. Re-call `parse()` for each pass.
- Use `read()` only when the file holds exactly one record; it raises loudly on zero or multiple records.
- Duplicate IDs make `to_dict`/`index`/`index_db` raise `ValueError: Duplicate key`. Provide a `key_function` to derive unique keys.
- `record.description` is the whole header including the id; `record.id` is just its first whitespace token.
- Mention if your file is compressed (.gz) so the agent opens it in text mode; only BGZF-compressed files support indexed random access.
- For FASTQ, tell the agent the quality encoding when known; the wrong variant string can silently shift every score by 31.

## Related Skills

- write-sequences - Write parsed sequences to new files
- filter-sequences - Filter sequences by criteria after reading
- format-conversion - Convert between formats
- compressed-files - Read gzip/bzip2/BGZF compressed sequence files
- fastq-quality - FASTQ quality encoding and offset selection
- sequence-manipulation/seq-objects - Work with parsed SeqRecord objects
- database-access/entrez-fetch - Fetch sequences from NCBI instead of local files

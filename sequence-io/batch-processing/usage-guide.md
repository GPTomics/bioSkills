# Batch Processing - Usage Guide

## Overview

This skill enables AI agents to help you process many sequence files at once (count, merge, split, convert, summarize) without exhausting RAM, choosing between Biopython streaming, on-disk indexing, and faster readers (pysam, pyfastx) based on data size.

## Prerequisites

```bash
pip install biopython
pip install pysam pyfastx   # optional, faster readers for huge FASTQ
```

## Quick Start

Tell your AI agent what you want to do:

- "Merge all FASTA files in the data folder into one file"
- "Split this large file into 1000-record chunks"
- "Count sequences in each file in the directory without loading them"
- "Build a random-access index across all my FASTA files"
- "Convert all GenBank files to FASTA"
- "Summarize counts and lengths for every file as a CSV"

## Example Prompts

### Counting
> "Count reads in every .fastq file under reads/ without loading them into memory"

### Merging
> "Combine all .fasta files in samples/ into one file, tagging each record with its source filename"

### Splitting
> "Split large.fasta into files of 500 sequences each"

### Random Access
> "Build a persistent index across all FASTA files in data/ so I can look up records by id later"

### Batch Conversion
> "Convert all .gb files in genbank/ to FASTA in fasta/"

### Statistics
> "Generate a summary CSV with sequence counts and min/max/avg lengths for all FASTA files"

### Large Data
> "Iterate over a 50-million-read FASTQ as fast as possible to tally read lengths"

## What the Agent Will Do
1. Choose a reader by data size: SeqIO for rich records, index_db for random access, pysam/pyfastx for huge FASTQ
2. Stream records with generators so peak memory stays at one record, never the whole directory
3. Apply the requested operation (count, merge, split, convert, summarize) per file
4. Write results and clean up index/output artifacts

## Tips

- Never `list(SeqIO.parse(...))` a large file; iterate the generator or use `sum(1 for _ in ...)` to count
- `SeqIO.parse` is one-pass and exhausts silently; re-create it per pass or use `index_db` for repeated lookups
- Use `SeqIO.index_db()` (not `to_dict()`) for random access across many or huge files; the SQLite index persists across sessions
- Duplicate ids across merged files raise `ValueError: Duplicate key` in `index_db`/`to_dict`; prefix ids by filename to keep them unique
- `SeqIO` is slow over tens of millions of reads because it builds a full record each; switch to `pysam.FastxFile` (linear) or `pyfastx` (indexed)
- Only BGZF-compressed files are indexable; plain gzip is not seekable, so re-`bgzip` before indexing
- Use `multiprocessing.Pool` for CPU-bound per-file work, `ThreadPoolExecutor` for I/O-bound work
- Use `rglob()` for recursive directory search and `pathlib.Path` for cross-platform paths

## Related Skills

- read-sequences - parse, index, and index_db semantics for each file
- filter-sequences - apply per-record filters while streaming a batch
- sequence-statistics - N50 and length distributions across files
- format-conversion - batch format conversion and its data-loss traps
- compressed-files - BGZF vs plain gzip for indexable batch random access
- paired-end-fastq - keep R1/R2 synchronized when batch-filtering mates
- database-access/entrez-fetch - batch download sequences from NCBI

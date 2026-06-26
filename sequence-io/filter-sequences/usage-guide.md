# Filter Sequences - Usage Guide

## Overview

This skill enables AI agents to help you filter and select sequences by length, ID, GC content, N content, motifs, regex patterns, and description using Biopython. Filtering streams record by record so large FASTA/FASTQ files never load into memory. For paired R1/R2 reads, the agent routes to the paired-end-fastq skill so mates stay synchronized.

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:

- "Keep only sequences longer than 500 bp"
- "Extract sequences matching these IDs"
- "Filter out sequences with N's"
- "Select sequences with GC content between 40 and 60 percent"
- "Drop reads shorter than 100 bp and with more than 5% N's"

## Example Prompts

### By Length
> "Remove sequences shorter than 200 bp from my FASTA file"

### By ID
> "Extract sequences with IDs listed in wanted_ids.txt and exclude everything else"

### By Content
> "Filter out sequences containing more than 5% N bases"

### By GC Content
> "Keep only sequences with GC content above 50%, counting N's against the full length"

### By Motif or Pattern
> "Keep sequences that contain an EcoRI site, ignoring soft-masked lowercase"

### Combined
> "Filter sequences: length >= 100, no N's, GC between 40 and 60 percent"

### Paired-end
> "Filter my R1/R2 FASTQ pair by quality and keep the mates synchronized"

## What the Agent Will Do
1. Check whether the input is single-file or paired R1/R2 (paired routes to paired-end-fastq).
2. Stream records from the input file one at a time.
3. Apply the requested criteria, setting the GC `ambiguous=` mode and uppercasing for content matches.
4. Write survivors to the output file and report the count.

## Tips

- Use generator expressions for large files; only random sampling and splitting need to load all records.
- For ID-based filtering, load IDs into a set for O(1) lookup.
- `gc_fraction()` returns a fraction 0-1, not a percent - use 0.4, not 40, for thresholds.
- Set `gc_fraction(seq, ambiguous=...)` explicitly: `'remove'` drops N's, `'ignore'` dilutes GC, `'weighted'` adds each code's expected GC.
- `Seq` preserves case - `.upper()` before motif/regex/uppercase tests so soft-masked lowercase bases are not missed.
- `SeqIO.parse()` is one-pass; recreate the generator if you need a second loop.
- Never filter one mate of a paired-end set alone - it desyncs R1/R2 and silently mismaps downstream.

## Related Skills

- read-sequences - Parse sequences before filtering
- write-sequences - Write filtered sequences to output
- fastq-quality - Filter FASTQ by per-base quality scores and encoding
- paired-end-fastq - Synchronized filtering of R1/R2 with orphan handling
- sequence-manipulation/sequence-properties - Per-sequence GC, length, and composition
- sequence-manipulation/motif-search - Filter by complex motif patterns
- alignment-files/alignment-filtering - Filter aligned reads with samtools view -f/-F

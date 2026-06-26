# Compressed Files - Usage Guide

## Overview

This skill enables AI agents to help you read, write, and randomly access compressed sequence files (.gz, .bz2, .xz, .bgz) using Biopython plus bgzip/samtools. The central idea: any of these formats can be read sequentially, but only BGZF (blocked gzip) is seekable, so only BGZF can be indexed for random access by `SeqIO.index`, `samtools faidx`, or `tabix`.

## Prerequisites

```bash
pip install biopython
# Optional, for CLI conversion and faidx indexing:
#   conda install -c bioconda htslib samtools   (provides bgzip + samtools)
```

gzip, bz2, and lzma are built into Python; no install needed for sequential read/write.

## Quick Start

Tell your AI agent what you want to do:

- "Read my gzipped FASTQ file"
- "Save these sequences to a gzip-compressed FASTA"
- "Convert reads.fastq.gz to an indexable BGZF file"
- "Randomly pull gene_042 out of a large compressed FASTA"
- "Why does samtools faidx reject my .gz file?"
- "Count sequences in reads.fastq.gz without decompressing to disk"

## Example Prompts

### Reading
> "Parse the gzipped FASTA and print each sequence id and length."
> "My .bz2 FASTQ won't parse with SeqIO - what handle mode do I need?"

### Writing
> "Write these records straight to a gzip-compressed FASTA without an intermediate plain file."
> "Save this FASTA as BGZF so I can index it later."

### Random access and indexing
> "Convert sequences.fasta.gz to BGZF and extract the region gene_042:1-200."
> "Build a persistent SQLite index over several BGZF FASTA files and fetch records by id."
> "samtools faidx says my file is 'not compressed with bgzip' - how do I fix it?"

### Conversion
> "Re-compress this plain gzip genome as BGZF without writing a decompressed copy to disk."
> "Decompress reads.fastq.gz back to plain reads.fastq."

## What the Agent Will Do
1. Identify the real compression (suffix is not authoritative - a `.gz` may be plain gzip or BGZF; `bgzip -t` confirms BGZF).
2. Open the file with the matching handler in TEXT mode (`gzip.open`/`bz2.open`/`lzma.open`/`bgzf.open`, `'rt'` or `'wt'`).
3. For random access, re-compress to BGZF (`zcat | bgzip` or the Bio.bgzf writer) and index with `SeqIO.index`/`index_db` or `samtools faidx`.
4. Stream records rather than loading them all into RAM for large files.

## Tips

- Always use TEXT mode (`'rt'`/`'wt'`) with SeqIO; a binary `'rb'` handle raises `TypeError`.
- A plain `.gz` is NOT seekable - `SeqIO.index` and faidx require BGZF; convert with `bgzip` first.
- `bgzf.open` reads BGZF input only; use `gzip.open` for plain gzip.
- A BGZF file is still valid gzip, so `gunzip`/`zcat` read it fine - the asymmetry only bites indexing tools.
- `samtools faidx` on a BGZF FASTA writes both `.fai` and `.gzi`; keep both or region extraction breaks.
- bzip2 and xz compress tighter than gzip but are slower and never seekable - use them for cold archives, not for files you index.
- For random access into a plain gzip you cannot re-bgzip, reach for `pyfastx` (seek-point index over the gzip stream).
- Virtual offsets from `bgzf.tell()` can be compared and re-seeked but never subtracted to get a byte length.

## Related Skills

- read-sequences - parse vs index vs index_db trade-offs for compressed handles
- write-sequences - write records through a compression handle
- batch-processing - stream many compressed files without loading them into RAM
- filter-sequences - keep paired reads in sync when filtering gzipped FASTQ
- alignment-files/sam-bam-basics - BAM is BGZF natively; samtools manages the compression

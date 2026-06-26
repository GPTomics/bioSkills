---
name: bio-batch-processing
description: Process many sequence files in batch (count, merge, split, convert, summarize) with memory-safe streaming and on-disk indexing using Biopython, pysam, or pyfastx. Use when iterating over a directory of FASTA/FASTQ files, merging or splitting datasets, building random access across many or huge files, or automating per-file operations without exhausting RAM.
tool_type: python
primary_tool: Bio.SeqIO
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+ (alternatives: pysam 0.22+, pyfastx 2.0+)

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Batch Processing

**"Process all my sequence files in a directory"** -> Iterate, merge, split, convert, and summarize across multiple sequence files without loading everything into RAM.
- Python: `SeqIO.parse()` + `Path.glob()` (BioPython, pathlib) for streaming
- Python: `SeqIO.index_db()` (BioPython) for persistent random access across many files
- Python: `pysam.FastxFile` (pysam) or `pyfastx` for fast iteration over huge FASTQ

## The Governing Principle

`list(SeqIO.parse(...))` materializes every `SeqRecord` in RAM at once. On a directory of large files this causes OOM. `SeqIO.parse()` itself returns a generator that holds one record at a time, so streaming is the default for batch work: iterate, never `list()`, unless the file is known-small and needs multiple passes.

For random access across many or huge files, do not load them. `SeqIO.index_db()` builds one on-disk SQLite index over a list of files that persists across sessions. That, not `to_dict()`, is the batch random-access tool.

For tens of millions of reads, `SeqIO` is slow by design: it constructs a full `SeqRecord` (a `Seq`, id/name/description, and a `letter_annotations` dict of per-base qualities) for every read. When the job is plain linear iteration, a thinner reader wins.

## Choosing a Reader

| Reader | Per-record object | Random access | Best for |
|--------|-------------------|---------------|----------|
| `Bio.SeqIO.parse` | full `SeqRecord` (rich API) | no (one-pass generator) | small/medium data needing the Biopython record API |
| `Bio.SeqIO.index_db` | reparsed `SeqRecord` on access | yes, on-disk SQLite, multi-file, persists | batch random access across many/huge files |
| `pysam.FastxFile` | thin entry (`.name/.sequence/.comment/.quality`) | no (linear, gzip sequential) | fast linear iteration over huge FASTQ |
| `pyfastx` | tuple/object via SQLite index | yes, into plain or gzipped FASTA/Q | random access + indexed reuse of gzipped files |

`pysam.FastxFile` exposes `.name`, `.sequence`, `.comment`, `.quality`, and `.get_quality_array()` (offset-removed int Phred, but it always subtracts 33, so it is correct only for Phred+33 data - for legacy Phred+64/Solexa stay on `SeqIO` with the explicit variant string). `pyfastx` builds a persistent `.fxi`/`.fqi` SQLite index and reads random records out of plain or gzipped files without re-bgzipping.

## Required Imports

```python
from pathlib import Path
from Bio import SeqIO
```

## Iterate and Count Across Files

Count by iterating, never by building a list. `len(list(SeqIO.parse(f)))` loads the whole file; `sum(1 for _ in ...)` holds one record at a time.

```python
for fasta_file in Path('data/').glob('*.fasta'):
    count = sum(1 for _ in SeqIO.parse(fasta_file, 'fasta'))
    print(f'{fasta_file.name}: {count} sequences')
```

Recursive search uses `rglob`:

```python
for gb_file in Path('data/').rglob('*.gb'):
    print(f'Found: {gb_file}')
```

For huge FASTQ where only sequence content matters, skip `SeqRecord` construction entirely:

```python
import pysam

with pysam.FastxFile('reads.fastq.gz') as fh:
    count = sum(1 for _ in fh)
```

## Random Access Across Many Files

**Goal:** Look up records by id across a whole directory of files, repeatedly, without holding them in RAM.

**Approach:** Build one persistent on-disk SQLite index over the file list with `index_db`. Reopen later with just the index path; lookups reparse single records from disk on demand.

**Reference (BioPython 1.83+):**

```python
from pathlib import Path
from Bio import SeqIO

files = [str(p) for p in Path('data/').glob('*.fasta')]
records = SeqIO.index_db('combined.idx', files, 'fasta')

print(len(records))            # total across all files
record = records['seq_00042']  # random access by id
records.close()
```

The index file persists. A later session calls `SeqIO.index_db('combined.idx')` with no file list and reopens instantly. Ids must be unique across the merged set: a collision raises `ValueError: Duplicate key`. `index_db` also indexes BGZF-compressed files; plain gzip is not seekable and cannot be indexed.

## Merge Files

**Goal:** Concatenate sequences from many files into one output without loading them all.

**Approach:** Chain per-file generators with `yield from` and stream straight into `SeqIO.write`, which consumes the generator one record at a time.

**Reference (BioPython 1.83+):**

```python
def all_records(directory, pattern, format):
    for filepath in Path(directory).glob(pattern):
        yield from SeqIO.parse(filepath, format)

count = SeqIO.write(all_records('data/', '*.fasta', 'fasta'), 'merged.fasta', 'fasta')
print(f'Merged {count} records')
```

### Merge with Source Tracking

**Goal:** Combine sequences from multiple files, tagging each record with its source filename.

**Approach:** Stream records through a generator that appends source metadata to the description before writing.

**Reference (BioPython 1.83+):**

```python
def records_with_source(directory, pattern, format):
    for filepath in Path(directory).glob(pattern):
        for record in SeqIO.parse(filepath, format):
            record.description = f'{record.description} [source={filepath.name}]'
            yield record

SeqIO.write(records_with_source('data/', '*.fasta', 'fasta'), 'merged_tracked.fasta', 'fasta')
```

When merging files that may share ids, decide upfront: write-then-merge tolerates duplicates (FASTA allows repeated ids), but any later `index_db`/`to_dict` over the merged file raises on the duplicate.

## Split Files

### Split by Number of Records

**Goal:** Divide a large file into chunks of N records each.

**Approach:** Consume the parse generator in fixed-size batches with `islice`, writing each batch to a numbered file. `islice` pulls only N records into memory per chunk, so an arbitrarily large input streams safely.

**Reference (BioPython 1.83+):**

```python
from itertools import islice

def split_file(input_file, format, records_per_file, output_prefix):
    records = SeqIO.parse(input_file, format)
    file_num = 1
    while True:
        batch = list(islice(records, records_per_file))
        if not batch:
            break
        output_file = f'{output_prefix}_{file_num}.{format}'
        SeqIO.write(batch, output_file, format)
        print(f'Wrote {len(batch)} records to {output_file}')
        file_num += 1

split_file('large.fasta', 'fasta', 1000, 'split')
```

On Python 3.12+, `itertools.batched(records, records_per_file)` yields the same fixed-size tuples without the manual `while`/`islice` loop.

### Split by Sequence ID Prefix

**Goal:** Group sequences into separate files by a shared id prefix (sample or chromosome).

**Approach:** Route each record to a per-prefix open output handle while streaming, so no group is fully held in RAM.

**Reference (BioPython 1.83+):**

```python
handles = {}
for record in SeqIO.parse('input.fasta', 'fasta'):
    prefix = record.id.split('_')[0]
    if prefix not in handles:
        handles[prefix] = open(f'{prefix}.fasta', 'w')
    SeqIO.write(record, handles[prefix], 'fasta')

for handle in handles.values():
    handle.close()
```

## Batch Convert

```python
for gb_file in Path('genbank/').glob('*.gb'):
    fasta_file = Path('fasta/') / gb_file.with_suffix('.fasta').name
    count = SeqIO.convert(str(gb_file), 'genbank', str(fasta_file), 'fasta')
    print(f'{gb_file.name} -> {fasta_file.name}: {count} records')
```

`SeqIO.convert` streams internally and never loads the whole file. GenBank-to-FASTA silently drops features, annotations, and qualifiers (FASTA stores only id, description, and sequence); see sequence-io/format-conversion before converting away annotated formats.

## Parallel Processing

For CPU-bound per-file work, distribute whole files across processes. Each worker streams its own file, so peak memory is one file's records per process, not the whole directory.

```python
from multiprocessing import Pool

def process_file(filepath):
    total = 0
    bp = 0
    for record in SeqIO.parse(filepath, 'fasta'):
        total += 1
        bp += len(record.seq)
    return {'file': filepath.name, 'count': total, 'total_bp': bp}

files = list(Path('data/').glob('*.fasta'))
with Pool(4) as pool:
    results = pool.map(process_file, files)
```

Use `concurrent.futures.ThreadPoolExecutor` instead for I/O-bound work (gzip decode, network filesystems); the GIL makes threads pointless for CPU-bound parsing.

## Summary Statistics

**Goal:** Build a per-file CSV of counts and length stats for a directory.

**Approach:** Stream each file once, accumulating count, total, min, and max as integers rather than collecting a length list per file.

**Reference (BioPython 1.83+):**

```python
import csv

summaries = []
for fasta_file in Path('data/').glob('*.fasta'):
    count = total = 0
    min_len = None
    max_len = 0
    for record in SeqIO.parse(fasta_file, 'fasta'):
        n = len(record.seq)
        count += 1
        total += n
        max_len = max(max_len, n)
        min_len = n if min_len is None else min(min_len, n)
    summaries.append({'file': fasta_file.name, 'sequences': count, 'total_bp': total,
                      'min_len': min_len or 0, 'max_len': max_len,
                      'avg_len': total / count if count else 0})

with open('summary.csv', 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=summaries[0].keys())
    writer.writeheader()
    writer.writerows(summaries)
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| `MemoryError` / process killed on a directory | `list(SeqIO.parse(...))` materializes every record at once | Stream the generator; iterate or `sum(1 for _ in ...)`; never `list()` a large file |
| Counting/merge job runs for minutes on tens of millions of reads | `SeqIO` builds a full `SeqRecord` per read | Use `pysam.FastxFile` for linear iteration, or `pyfastx` for indexed access |
| `ValueError: Duplicate key` from `index_db`/`to_dict` | Same id appears in more than one merged file | Make ids unique (prefix by filename) or supply a `key_function` |
| Second loop over `SeqIO.parse(...)` yields nothing | The generator is one-pass and exhausts silently | Re-create the generator per pass, or use `index_db` for repeated access |
| `index_db` fails on a `.gz` file | Plain gzip is not seekable | Re-compress with `bgzip`; only BGZF is indexable (sequence-io/compressed-files) |
| Annotations missing after batch convert | GenBank-to-FASTA drops all features silently | Keep an annotated format, or extract needed qualifiers first |

## Related Skills

- read-sequences - parse, index, and index_db semantics for each file
- filter-sequences - apply per-record filters while streaming a batch
- sequence-statistics - N50 and length distributions across files
- format-conversion - batch format conversion and its data-loss traps
- compressed-files - BGZF vs plain gzip for indexable batch random access
- paired-end-fastq - keep R1/R2 synchronized when batch-filtering mates
- database-access/entrez-fetch - batch download sequences from NCBI

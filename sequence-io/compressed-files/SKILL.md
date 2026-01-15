---
name: bio-compressed-files
description: Read and write compressed sequence files (gzip, bzip2) using Biopython with Python's compression modules. Use when working with .gz or .bz2 sequence files.
---

# Compressed Files

Handle gzip and bzip2 compressed sequence files with Biopython.

## Required Imports

```python
import gzip
import bz2
from Bio import SeqIO
```

## Reading Compressed Files

### Gzip (.gz)
```python
with gzip.open('sequences.fasta.gz', 'rt') as handle:
    for record in SeqIO.parse(handle, 'fasta'):
        print(record.id, len(record.seq))
```

**Important:** Use `'rt'` (read text) mode, not `'rb'` (read binary).

### Bzip2 (.bz2)
```python
with bz2.open('sequences.fasta.bz2', 'rt') as handle:
    for record in SeqIO.parse(handle, 'fasta'):
        print(record.id, len(record.seq))
```

## Writing Compressed Files

### Gzip (.gz)
```python
with gzip.open('output.fasta.gz', 'wt') as handle:
    SeqIO.write(records, handle, 'fasta')
```

### Bzip2 (.bz2)
```python
with bz2.open('output.fasta.bz2', 'wt') as handle:
    SeqIO.write(records, handle, 'fasta')
```

## Code Patterns

### Read Gzipped FASTQ
```python
with gzip.open('reads.fastq.gz', 'rt') as handle:
    records = list(SeqIO.parse(handle, 'fastq'))
print(f'Loaded {len(records)} reads')
```

### Count Records in Gzipped File
```python
with gzip.open('sequences.fasta.gz', 'rt') as handle:
    count = sum(1 for _ in SeqIO.parse(handle, 'fasta'))
print(f'{count} sequences')
```

### Convert Compressed to Uncompressed
```python
with gzip.open('input.fasta.gz', 'rt') as in_handle:
    records = SeqIO.parse(in_handle, 'fasta')
    SeqIO.write(records, 'output.fasta', 'fasta')
```

### Convert Uncompressed to Compressed
```python
records = SeqIO.parse('input.fasta', 'fasta')
with gzip.open('output.fasta.gz', 'wt') as out_handle:
    SeqIO.write(records, out_handle, 'fasta')
```

### Process Large Gzipped File (Memory Efficient)
```python
with gzip.open('large.fastq.gz', 'rt') as handle:
    for record in SeqIO.parse(handle, 'fastq'):
        if len(record.seq) >= 100:
            process(record)
```

### Auto-Detect Compression
```python
from pathlib import Path

def open_sequence_file(filepath, format):
    filepath = Path(filepath)
    if filepath.suffix == '.gz':
        handle = gzip.open(filepath, 'rt')
    elif filepath.suffix == '.bz2':
        handle = bz2.open(filepath, 'rt')
    else:
        handle = open(filepath, 'r')
    return SeqIO.parse(handle, format)
```

### Compress Existing File
```python
import shutil

with open('sequences.fasta', 'rb') as f_in:
    with gzip.open('sequences.fasta.gz', 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
```

## File Extensions

| Extension | Compression | Module |
|-----------|-------------|--------|
| `.gz` | gzip | `gzip` |
| `.bz2` | bzip2 | `bz2` |
| `.xz` | LZMA | `lzma` |

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `TypeError: a bytes-like object is required` | Used 'rb' mode | Use 'rt' for text mode |
| `UnicodeDecodeError` | Wrong encoding | Try `gzip.open(file, 'rt', encoding='latin-1')` |
| `gzip.BadGzipFile` | Not a gzip file | Check file extension matches actual format |
| `OSError: Not a gzipped file` | Corrupt or wrong format | Verify file integrity |

## Performance Notes

- Gzipped files cannot be indexed with SeqIO.index()
- For random access to compressed files, consider uncompressing first
- Bzip2 has better compression but slower read/write than gzip
- For very large files, process as stream rather than loading all into memory

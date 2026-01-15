---
name: bio-read-sequences
description: Read biological sequence files (FASTA, FASTQ, GenBank, EMBL) using Biopython Bio.SeqIO. Use when parsing sequence files, iterating multi-sequence files, or loading sequences for random access.
---

# Read Sequences

Read biological sequence data from files using Biopython's Bio.SeqIO module.

## Required Import

```python
from Bio import SeqIO
```

## Core Functions

### SeqIO.parse() - Multiple Records (Iterator)
Use for files with one or more sequences. Returns an iterator of SeqRecord objects.

```python
for record in SeqIO.parse('sequences.fasta', 'fasta'):
    print(record.id, len(record.seq))
```

**Important:** Always specify the format explicitly as the second argument.

### SeqIO.read() - Single Record Only
Use when file contains exactly one sequence. Raises error if zero or multiple records.

```python
record = SeqIO.read('single.fasta', 'fasta')
```

### SeqIO.to_dict() - Load All Into Memory
Use for random access by record ID. Loads entire file into memory.

```python
records = SeqIO.to_dict(SeqIO.parse('sequences.fasta', 'fasta'))
seq = records['sequence_id'].seq
```

### SeqIO.index() - Large File Random Access
Use for large files when you need random access without loading everything into memory.

```python
records = SeqIO.index('large.fasta', 'fasta')
seq = records['sequence_id'].seq
records.close()
```

## Common Formats

| Format | String | Typical Extension |
|--------|--------|-------------------|
| FASTA | `'fasta'` | .fasta, .fa, .fna, .faa |
| FASTQ | `'fastq'` | .fastq, .fq |
| GenBank | `'genbank'` or `'gb'` | .gb, .gbk |
| EMBL | `'embl'` | .embl |
| Swiss-Prot | `'swiss'` | .dat |
| PHYLIP | `'phylip'` | .phy |
| Clustal | `'clustal'` | .aln |
| Stockholm | `'stockholm'` | .sto |

## SeqRecord Object Attributes

After parsing, each record has these key attributes:

```python
record.id          # Sequence identifier (string)
record.name        # Sequence name (string)
record.description # Full description line (string)
record.seq         # Sequence data (Seq object)
record.features    # List of SeqFeature objects (GenBank/EMBL)
record.annotations # Dictionary of annotations
record.letter_annotations  # Per-letter annotations (quality scores)
record.dbxrefs     # Database cross-references
```

## Code Patterns

### Collect All Sequences Into a List
```python
records = list(SeqIO.parse('sequences.fasta', 'fasta'))
```

### Count Records Without Loading All
```python
count = sum(1 for _ in SeqIO.parse('sequences.fasta', 'fasta'))
```

### Get Sequence IDs Only
```python
ids = [record.id for record in SeqIO.parse('sequences.fasta', 'fasta')]
```

### Extract Sequences as Strings
```python
sequences = [str(record.seq) for record in SeqIO.parse('sequences.fasta', 'fasta')]
```

### Read GenBank with Features
```python
for record in SeqIO.parse('sequence.gb', 'genbank'):
    for feature in record.features:
        if feature.type == 'CDS':
            print(feature.qualifiers.get('product', ['Unknown'])[0])
```

### Access FASTQ Quality Scores
```python
for record in SeqIO.parse('reads.fastq', 'fastq'):
    qualities = record.letter_annotations['phred_quality']
    avg_quality = sum(qualities) / len(qualities)
```

### Read From File Handle
```python
with open('sequences.fasta', 'r') as handle:
    for record in SeqIO.parse(handle, 'fasta'):
        print(record.id)
```

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `ValueError: More than one record` | Used `read()` on multi-record file | Use `parse()` instead |
| `ValueError: No records found` | Used `read()` on empty file | Check file exists and has content |
| `ValueError: unknown format` | Typo in format string | Check format string spelling |
| `UnicodeDecodeError` | Binary file or wrong encoding | Open with `encoding='latin-1'` or check file |

## Decision Tree

```
Need to read sequences?
├── Single record in file?
│   └── Use SeqIO.read()
├── Multiple records?
│   ├── Need all in memory at once?
│   │   └── Use list(SeqIO.parse()) or SeqIO.to_dict()
│   ├── Process one at a time (memory efficient)?
│   │   └── Use SeqIO.parse() iterator
│   └── Large file, need random access by ID?
│       └── Use SeqIO.index()
```

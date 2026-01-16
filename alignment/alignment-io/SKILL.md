---
name: bio-alignment-io
description: Read, write, and convert multiple sequence alignment files using Biopython Bio.AlignIO. Supports Clustal, PHYLIP, Stockholm, FASTA, Nexus, and other alignment formats for phylogenetics and conservation analysis.
tool_type: python
primary_tool: Bio.AlignIO
---

# Alignment File I/O

Read, write, and convert multiple sequence alignment files in various formats.

## Required Import

```python
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
```

## Supported Formats

| Format | Extension | Read | Write | Description |
|--------|-----------|------|-------|-------------|
| `clustal` | .aln | Yes | Yes | Clustal W/X output |
| `fasta` | .fasta, .fa | Yes | Yes | Aligned FASTA |
| `phylip` | .phy | Yes | Yes | Interleaved PHYLIP |
| `phylip-sequential` | .phy | Yes | Yes | Sequential PHYLIP |
| `phylip-relaxed` | .phy | Yes | Yes | PHYLIP with long names |
| `stockholm` | .sto, .stk | Yes | Yes | Pfam/Rfam annotated |
| `nexus` | .nex | Yes | Yes | NEXUS format |
| `emboss` | .txt | Yes | No | EMBOSS tools output |
| `fasta-m10` | .txt | Yes | No | FASTA -m 10 output |
| `maf` | .maf | Yes | Yes | Multiple Alignment Format |
| `mauve` | .xmfa | Yes | No | progressiveMauve output |
| `msf` | .msf | Yes | No | GCG MSF format |

## Reading Alignments

### Single Alignment File
```python
from Bio import AlignIO

alignment = AlignIO.read('alignment.aln', 'clustal')
print(f'Alignment length: {alignment.get_alignment_length()}')
print(f'Number of sequences: {len(alignment)}')
```

### Multiple Alignments in One File
```python
for alignment in AlignIO.parse('multi_alignment.sto', 'stockholm'):
    print(f'Alignment with {len(alignment)} sequences, length {alignment.get_alignment_length()}')
```

### Read as List
```python
alignments = list(AlignIO.parse('alignments.phy', 'phylip'))
print(f'Read {len(alignments)} alignments')
```

## Writing Alignments

### Write Single Alignment
```python
AlignIO.write(alignment, 'output.fasta', 'fasta')
```

### Write Multiple Alignments
```python
alignments = [alignment1, alignment2, alignment3]
count = AlignIO.write(alignments, 'output.sto', 'stockholm')
print(f'Wrote {count} alignments')
```

### Write to Handle
```python
with open('output.aln', 'w') as handle:
    AlignIO.write(alignment, handle, 'clustal')
```

## Format Conversion

### Direct Conversion (Most Efficient)
```python
AlignIO.convert('input.aln', 'clustal', 'output.phy', 'phylip')
```

### With Alphabet Specification
```python
AlignIO.convert('input.sto', 'stockholm', 'output.nex', 'nexus', molecule_type='DNA')
```

### Manual Conversion (When Modification Needed)
```python
alignment = AlignIO.read('input.aln', 'clustal')
# ... modify alignment ...
AlignIO.write(alignment, 'output.fasta', 'fasta')
```

## Accessing Alignment Data

```python
alignment = AlignIO.read('alignment.aln', 'clustal')

# Iterate over sequences
for record in alignment:
    print(f'{record.id}: {record.seq}')

# Access by index
first_seq = alignment[0]
last_seq = alignment[-1]

# Slice columns
column_slice = alignment[:, 10:20]  # Columns 10-19

# Get specific column
column = alignment[:, 5]  # Column 5 as string
```

## Working with Alignment Objects

### Get Alignment Properties
```python
alignment = AlignIO.read('alignment.aln', 'clustal')

length = alignment.get_alignment_length()
num_seqs = len(alignment)
seq_ids = [record.id for record in alignment]
```

### Slice Alignments
```python
# Get subset of sequences
subset = alignment[0:5]  # First 5 sequences

# Get subset of columns
trimmed = alignment[:, 50:150]  # Columns 50-149

# Combine slicing
region = alignment[0:5, 50:150]  # 5 sequences, columns 50-149
```

## Creating Alignments Programmatically

```python
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

records = [
    SeqRecord(Seq('ACTGACTGACTG'), id='seq1'),
    SeqRecord(Seq('ACTGACT-ACTG'), id='seq2'),
    SeqRecord(Seq('ACTG-CTGACTG'), id='seq3'),
]
alignment = MultipleSeqAlignment(records)
AlignIO.write(alignment, 'new_alignment.fasta', 'fasta')
```

## Format-Specific Notes

### PHYLIP Format
```python
# Standard PHYLIP (10 char names, interleaved)
alignment = AlignIO.read('file.phy', 'phylip')

# Sequential PHYLIP
alignment = AlignIO.read('file.phy', 'phylip-sequential')

# Relaxed PHYLIP (allows longer names)
alignment = AlignIO.read('file.phy', 'phylip-relaxed')
```

### Stockholm Format (with Annotations)
```python
alignment = AlignIO.read('pfam.sto', 'stockholm')

# Access annotations
for record in alignment:
    print(record.id, record.annotations)
```

### Clustal Format
```python
# Clustal preserves conservation symbols in file but not when parsed
alignment = AlignIO.read('clustal.aln', 'clustal')
```

## Batch Processing Multiple Files

```python
from pathlib import Path

input_dir = Path('alignments/')
output_dir = Path('converted/')

for input_file in input_dir.glob('*.aln'):
    alignment = AlignIO.read(input_file, 'clustal')
    output_file = output_dir / f'{input_file.stem}.fasta'
    AlignIO.write(alignment, output_file, 'fasta')
```

## Alternative: Bio.Align Module I/O

The newer `Bio.Align` module provides its own I/O functions that return `Alignment` objects (instead of `MultipleSeqAlignment`). These support additional formats and provide access to modern alignment features.

```python
from Bio import Align

# Read single alignment (returns Alignment object)
alignment = Align.read('alignment.aln', 'clustal')

# Parse multiple alignments
for alignment in Align.parse('multi.sto', 'stockholm'):
    print(f'Alignment with {len(alignment)} sequences')

# Write alignment
Align.write(alignment, 'output.fasta', 'fasta')
```

### When to Use Which

| Use Case | Module |
|----------|--------|
| Legacy code, MultipleSeqAlignment needed | `Bio.AlignIO` |
| Modern features (counts, substitutions) | `Bio.Align` |
| Format conversion | Either works |
| Working with pairwise alignments | `Bio.Align` |

## Quick Reference: Common Operations

| Task | Code |
|------|------|
| Read single alignment | `AlignIO.read(file, format)` |
| Read multiple alignments | `AlignIO.parse(file, format)` |
| Write alignment(s) | `AlignIO.write(align, file, format)` |
| Convert format | `AlignIO.convert(in_file, in_fmt, out_file, out_fmt)` |
| Get length | `alignment.get_alignment_length()` |
| Get sequence count | `len(alignment)` |
| Slice columns | `alignment[:, start:end]` |

## Common Errors

| Error | Cause | Solution |
|-------|-------|----------|
| `ValueError: No records` | Empty file | Check file path and format |
| `ValueError: More than one record` | Multiple alignments with `read()` | Use `parse()` instead |
| `ValueError: Sequences different lengths` | Invalid alignment | Ensure all sequences same length |
| `ValueError: unknown format` | Unsupported format string | Check supported formats list |

## Related Skills

- **pairwise-alignment** - Create pairwise alignments with PairwiseAligner
- **msa-parsing** - Analyze alignment content and annotations
- **alignment-statistics** - Calculate conservation and identity
- **sequence-io/format-conversion** - Convert sequence (non-alignment) formats

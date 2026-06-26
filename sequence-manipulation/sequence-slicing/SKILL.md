---
name: bio-sequence-slicing
description: Slice, extract, and concatenate biological sequences and annotated records using Biopython. Use when extracting subsequences by position, splicing exons into a transcript, joining sequences, or carrying a sub-region of an annotated record (with quality scores and features) into a new record.
tool_type: python
primary_tool: Bio.Seq
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Sequence Slicing

Extract sub-regions, splice non-contiguous regions, and concatenate sequences and annotated records.

**"Extract a subsequence"** -> Slice a Seq with 0-based half-open coordinates.
- Python: `seq[start:end]` (Bio.Seq)

**"Pull out a sub-region but keep its quality scores and features"** -> Slice the SeqRecord, not the bare Seq.
- Python: `record[start:end]` (Bio.SeqRecord)

**"Splice exons into a transcript"** -> Extract each region and concatenate.
- Python: `sum((seq[s:e] for s, e in coords), Seq(''))` or the `+` operator

## The Governing Principle

Slicing a bare `Seq` is pure string math: `seq[start:end]` returns a new `Seq`, half-open, with no metadata to lose. Slicing a `SeqRecord` carries metadata, and the rule for WHAT survives is the single most error-prone part of this skill:

`record[start:end]` (verified against `Bio/SeqRecord.__getitem__`):
- PRESERVES `id`, `name`, `description`, and `molecule_type`.
- AUTO-SLICES `letter_annotations` (per-letter data such as PHRED `phred_quality`) to match the new coordinates -- this is why a FASTQ slice keeps the right per-base qualities for free.
- KEEPS only features FULLY CONTAINED in `[start:end]`; their locations are recalculated relative to the new start.
- SILENTLY DROPS the `annotations` dict (organism, taxonomy, references, comments), the `dbxrefs` list, and any feature that STRADDLES the slice boundary (dropped whole, never truncated). A non-trivial stride (`record[::2]`) drops features entirely.

Nothing warns when annotations vanish. The GenBank `source` feature spans the whole record, so it straddles almost any slice and disappears along with organism/taxonomy. To carry that metadata across, copy it explicitly:

```python
sub = record[start:end]
sub.annotations = record.annotations.copy()
```

and re-add any boundary-straddling feature manually (with a clamped, recalculated location) if a truncated copy is needed.

## Required Imports

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
```

## Coordinate Systems: the 0-based vs 1-based trap

Python and Biopython slicing is 0-based and half-open: `seq[start:end]` includes `start`, excludes `end`, and returns `end - start` letters. File formats disagree, and mixing them is a SILENT off-by-one (no error, just the wrong bases):

| Source | Convention | Position 1234..5678 means |
|--------|------------|---------------------------|
| Python / Bio.Seq slice | 0-based, half-open | `seq[1234:5678]` |
| GenBank / EMBL / GFF / VCF feature line | 1-based, INCLUSIVE | `seq[1233:5678]` (subtract 1 from start only) |
| BED file | 0-based, half-open | `seq[1234:5678]` (already matches Python) |

The asymmetry is the catch: convert a 1-based inclusive interval by subtracting 1 from the START only; the end already lands correctly because Python's exclusive end cancels the inclusive end. Reading a coordinate straight off a GFF and slicing `seq[start:end]` without the `-1` silently shifts everything one base left.

`Bio.SeqFeature` locations sidestep this entirely: they store a 0-based start and a Python-style end, so `int(feature.location.start):int(feature.location.end)` slices the parent directly, and `feature.extract(record.seq)` does the same automatically (handling strand and compound/joined locations).

```python
def extract_1based(seq, start, end):
    '''Extract a 1-based inclusive interval (GenBank/GFF style).'''
    return seq[start - 1:end]
```

## Slicing a Bare Seq

Slicing returns a `Seq` (not a string); negative indices and strides behave exactly like `str` (Seq has behaved like `str` since BioPython 1.78).

```python
seq = Seq('ATGCGATCGATCG')
seq[0]       # 'A'  single base, 0-indexed -> returns a str
seq[-1]      # 'G'  last base
seq[0:3]     # Seq('ATG')   first 3 bases
seq[-5:]     # Seq('GATCG')  last 5
seq[::2]     # Seq('AGGTGTG')  every 2nd base (stride)
seq[::-1]    # Seq('GCTAGCTAGCGTA')  reversed (not the reverse complement)
```

`str(record.seq)` returns the raw string, but raises `UndefinedSequenceError` when the record's sequence content is undefined (e.g. `Seq(None, length=n)` from a header-only FASTA or a pysam-backed record). Guard with `len()` (always defined) before forcing the content to a string.

## Code Patterns

### Splice Non-Contiguous Regions (Exons -> Transcript)

**Goal:** Join several separated regions of a genomic sequence into one continuous sequence.

**Approach:** Extract each region with half-open coordinates and concatenate. `sum()` needs an explicit `Seq('')` start value because the default `0` cannot be added to a `Seq`.

```python
def extract_regions(seq, regions):
    '''Concatenate multiple [start, end) regions in order.'''
    return sum((seq[start:end] for start, end in regions), Seq(''))

exon_coords = [(0, 50), (100, 150), (200, 250)]
mrna = extract_regions(genomic_seq, exon_coords)
```

For a real annotated transcript, let the feature do the work -- `feature.extract` honors strand and joined exon locations:

```python
for feature in record.features:
    if feature.type == 'mRNA':
        transcript = feature.extract(record.seq)
```

### Carry a Sub-Region into a New Annotated Record

**Goal:** Keep id, per-base quality, and contained features when extracting a window, and decide deliberately what metadata to carry.

**Approach:** Slice the `SeqRecord` (qualities and contained features ride along automatically), then explicitly copy the `annotations` dict, which slicing always drops.

```python
sub = record[100:400]                      # qualities + contained features auto-sliced
sub.annotations = record.annotations.copy()  # organism/taxonomy/refs would be lost otherwise
sub.id = f'{record.id}:101-400'            # 1-based label for humans
```

To build a fresh record from a bare Seq slice instead (no source metadata to carry):

```python
sub = SeqRecord(record.seq[100:400], id=f'{record.id}_sub', description='positions 101-400')
```

### Extract a Feature by Type

```python
for record in SeqIO.parse('sequence.gb', 'genbank'):
    for feature in record.features:
        if feature.type == 'CDS':
            cds = feature.extract(record.seq)      # strand-aware
            gene = feature.qualifiers.get('gene', ['?'])[0]
```

### Concatenate Sequences and Records

```python
seq1 + seq2                      # Seq + Seq -> Seq
seq1 + 'NNNN'                    # Seq + str -> Seq
Seq('NNN').join([s1, s2, s3])   # linker between each -> Seq
```

Adding `SeqRecord` objects works (`rec1 + rec2` concatenates sequences and per-letter annotations), but follows the same rule as slicing: the result keeps `id`/`name`/`description` only when both share them, and the `annotations` dict is reset. Set metadata on the result explicitly.

### Split into Codons or Fixed Chunks

```python
def split_codons(seq):
    '''Whole codons only; trailing 1-2 nt remainder is dropped.'''
    return [seq[i:i + 3] for i in range(0, len(seq) - len(seq) % 3, 3)]

def chunk_sequence(seq, size):
    '''Fixed-size chunks; final chunk may be shorter.'''
    return [seq[i:i + size] for i in range(0, len(seq), size)]
```

### Tile Overlapping Windows

```python
def sliding_windows(seq, window_size, step=1):
    for i in range(0, len(seq) - window_size + 1, step):
        yield i, seq[i:i + window_size]
```

### Flanking Region Around a Position

```python
def get_flanking(seq, position, flank):
    '''Clamp to sequence ends so the slice never runs past the edges.'''
    start = max(0, position - flank)
    end = min(len(seq), position + flank + 1)
    return seq[start:end]
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Organism/taxonomy/references gone from a sub-record | `record[start:end]` silently drops the `annotations` dict and `dbxrefs` | `sub.annotations = record.annotations.copy()` after slicing |
| A feature spanning the cut is missing from the slice | Features straddling the boundary are dropped whole, not truncated | Re-add manually with a clamped, recalculated location |
| All features gone after `record[::2]` | A non-trivial stride drops features entirely | Slice without a stride, or rebuild features by hand |
| Everything shifted one base left | GFF/GenBank 1-based start sliced as if 0-based | Subtract 1 from the START only: `seq[start-1:end]` |
| `UndefinedSequenceError` on `str(record.seq)` | Sequence content undefined (`Seq(None, length=n)`) | Use `len(record)`; do not force undefined content to a string |
| `TypeError` from `sum(slices)` | Default start `0` cannot add to a `Seq` | Pass a start: `sum(slices, Seq(''))` |
| Reversed but wrong strand | `seq[::-1]` reverses only; it does not complement | Use `seq.reverse_complement()` (see reverse-complement) |
| `IndexError` on single-base index | Position past the end | Check `len(seq)` first; slices clamp but `seq[i]` does not |

## Decision Guide

- Bare sequence, no metadata to keep -> slice the `Seq`: `seq[start:end]`.
- Need per-base quality or contained features to ride along -> slice the `SeqRecord`: `record[start:end]`, then copy `annotations`.
- Coordinates came from a GFF/GenBank/EMBL/VCF line -> subtract 1 from the start before slicing.
- Coordinates came from a BED file -> use as-is (already 0-based half-open).
- Strand-aware or joined/compound location -> `feature.extract(record.seq)`, never a manual slice.
- Joining separated regions -> `sum((seq[s:e] for s, e in coords), Seq(''))`.

## Related Skills

- seq-objects - Create Seq/SeqRecord objects and handle undefined sequence content
- reverse-complement - Reverse-complement an extracted region (slicing reverses but does not complement)
- transcription-translation - Translate an extracted CDS or spliced transcript
- sequence-io/read-sequences - Parse GenBank/FASTQ records (with features and qualities) to slice
- genome-intervals/gtf-gff-handling - Read 1-based GFF/GTF feature coordinates before slicing
- alignment-files/sam-bam-basics - Extract sequences from BAM regions with samtools

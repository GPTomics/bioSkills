---
name: bio-seq-objects
description: Create and manipulate Seq, MutableSeq, and SeqRecord objects using Biopython. Use when creating sequences from strings, modifying sequence data in-place, building annotated records for file output, or debugging post-1.78 Bio.Alphabet and immutability errors.
tool_type: python
primary_tool: Bio.Seq
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Seq Objects

Create and manipulate biological sequence objects using Biopython.

**"Create a sequence object"** -> Wrap a raw string in a typed sequence container for biological operations.
- Immutable: `Seq('ATGC')` (BioPython) - string-like, supports complement/translate
- Mutable: `MutableSeq('ATGC')` (BioPython) - supports in-place edits
- Annotated: `SeqRecord(Seq(...), id=...)` (BioPython) - adds metadata for file I/O

## The governing principle: no alphabet, no validation

`Bio.Alphabet` was removed entirely in Biopython 1.78 (2020-09-04). `Seq` and `SeqRecord` lost their `.alphabet` attribute, and any old-style construction fails LOUD: `from Bio.Alphabet import IUPAC` raises ImportError, and `Seq('ACGT', IUPAC.unambiguous_dna)` raises TypeError. Molecule type now lives as a SeqRecord annotation, `record.annotations['molecule_type'] = 'DNA'`, consumed by the GenBank/EMBL writers.

The consequence governs everything downstream: no `Seq` operation validates its alphabet anymore. A protein passed to `reverse_complement()` or `transcribe()` returns silent garbage rather than an error. Sibling skills (transcription-translation, reverse-complement) inherit this - the burden is on the caller to track what kind of molecule a `Seq` holds.

## Required Imports

```python
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord
```

## Core Objects

### Seq - Immutable Sequence

Immutable and behaves like `str` since 1.78: indexing, slicing, `+`, `*`, `.upper()`, `in`, `.count()`, `.find()` all work. In-place edits raise: `seq[0] = 'A'` -> TypeError (LOUD). Use MutableSeq for edits.

```python
seq = Seq('ATGCGATCGATCG')

len(seq)           # length
seq[0]             # first base
seq[0:10]          # slice (returns Seq)
str(seq)           # text form (see bytes note below)
'ATG' in seq       # membership test
seq.count('G')     # count occurrences
seq.find('ATG')    # position (-1 if not found)
seq.upper()        # uppercase (returns Seq)
seq * 3            # repeat
```

Since 1.79 `Seq` is backed by `bytes` (and `MutableSeq` by `bytearray`), NOT a `str` subclass. Use `str(seq)` for text and `bytes(seq)` for bytes. `isinstance(seq, str)` is always False - code that type-checks with `isinstance(x, str)` to detect sequences silently skips every `Seq`; test `isinstance(x, (Seq, MutableSeq))` instead.

### MutableSeq - Mutable Sequence

A `bytearray`-backed sequence for in-place editing; required when an operation needs `inplace=True`.

```python
mut_seq = MutableSeq('ATGCGATCG')
mut_seq[0] = 'C'              # modify single position
mut_seq[0:3] = 'GGG'          # replace slice
mut_seq.append('A')           # add to end
mut_seq.insert(0, 'G')        # insert at position
mut_seq.pop()                 # remove and return last
mut_seq.remove('G')           # remove first occurrence
mut_seq.reverse()             # reverse in place
```

Convert between types (a `MutableSeq` is unhashable and cannot be a dict key or used in `SeqIO.write`, so cast back to `Seq` when done editing):

```python
seq = Seq(mut_seq)            # MutableSeq -> Seq
mut_seq = MutableSeq(seq)     # Seq -> MutableSeq
```

### Undefined and partially-defined sequences

`UndefinedSequenceError` (added 1.79, a subclass of `ValueError`) models a sequence whose length is known but whose content is not - produced by lazy/partial file parsers. A `Seq(None, length=20)` reports `len() == 20` but raises on any attempt to read the bytes.

```python
undef = Seq(None, length=20)
len(undef)            # 20 - fine
str(undef)            # raises UndefinedSequenceError (subclass of ValueError)

partial = Seq({3: 'ACGT'}, length=10)   # only positions 3-6 defined
str(partial[3:7])     # 'ACGT' - defined region reads fine
str(partial)          # raises - undefined positions
```

Note: `complement()`/`reverse_complement()` on an undefined `Seq` return self rather than crash, but any read of the bytes raises. Guard reads of records from lazy parsers with `try`/`except UndefinedSequenceError` only where content access is genuinely optional.

### SeqRecord - Annotated Sequence

Sequence plus metadata for file I/O and analysis.

```python
record = SeqRecord(Seq('ATGCGATCG'), id='gene1', name='example_gene', description='An example gene sequence')

record.seq                 # the Seq object
record.id                  # identifier string
record.name                # name string
record.description         # description string
record.features            # list of SeqFeature objects
record.annotations         # dict (organism, molecule_type, topology, ...)
record.letter_annotations  # per-letter annotations (e.g. phred_quality)
record.dbxrefs             # database cross-references
```

### SeqRecord transformations

**Goal:** Transform whole records (reverse-complement, translate, slice) while keeping metadata coherent.

**Approach:** Use SeqRecord methods that return new records with features remapped to the new coordinate frame; pass `id`/`description` explicitly because they are NOT carried automatically.

```python
rc_record = record.reverse_complement(id=f'{record.id}_rc', description='reverse complement')
protein_record = record.translate(id=f'{record.id}_protein', to_stop=True)
fasta_str = record.format('fasta')      # quick in-memory file-format string
```

Slicing a SeqRecord remaps features but silently DROPS `annotations`, `dbxrefs`, and any feature that straddles a slice boundary - `subset = record[10:50]` returns a record with empty `annotations`. Re-attach `molecule_type` (and anything else a writer needs) on the slice before writing.

```python
subset = record[10:50]                          # features clipped; annotations dropped
subset.annotations['molecule_type'] = 'DNA'     # restore before GenBank/EMBL write
```

## Code Patterns

### Create Seq from String
```python
dna = Seq('ATGCGATCGATCG')
rna = Seq('AUGCGAUCGAUCG')
protein = Seq('MRCRS')
```

### Create SeqRecord with annotations for GenBank output
```python
record = SeqRecord(Seq('ATGCGATCG'), id='gene1', description='Example')
record.annotations['organism'] = 'Homo sapiens'
record.annotations['molecule_type'] = 'DNA'   # required by GenBank/EMBL writers
```

### Build SeqRecord with a feature
```python
from Bio.SeqFeature import SeqFeature, FeatureLocation

record = SeqRecord(Seq('ATGCGATCGATCG'), id='gene1')
feature = SeqFeature(FeatureLocation(0, 9), type='CDS', qualifiers={'product': ['Example protein']})
record.features.append(feature)
```

### Batch create SeqRecords
```python
sequences = ['ATGC', 'GCTA', 'TTAA']
records = [SeqRecord(Seq(s), id=f'seq_{i}') for i, s in enumerate(sequences)]
```

### Copy a SeqRecord
```python
from copy import deepcopy
new_record = deepcopy(record)   # deep copy; plain assignment shares features/annotations
new_record.id = 'modified_copy'
```

### Join sequences with a linker
```python
combined_seq = seq1 + Seq('NNNN') + seq2
combined_record = SeqRecord(combined_seq, id='combined')
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| `ImportError: No module named 'Bio.Alphabet'` (or `cannot import name 'IUPAC'`) | `Bio.Alphabet` removed in 1.78 | Drop the alphabet argument; set `record.annotations['molecule_type']` instead |
| `TypeError: 'Seq' object does not support item assignment` | Editing an immutable `Seq` in place | Use `MutableSeq`, or rebuild with slicing/concatenation |
| `UndefinedSequenceError` on `str(seq)`/`print(seq)` | Sequence from a lazy/partial parser (`Seq(None, length=n)`) has known length but no content | Avoid reading bytes, or guard with `except UndefinedSequenceError` (subclass of `ValueError`) |
| `isinstance(seq, str)` is False, type-check skips the sequence | Since 1.79 `Seq` is `bytes`-backed, not a `str` subclass | Test `isinstance(x, (Seq, MutableSeq))`; use `str(seq)` for text |
| `ValueError: missing molecule_type` writing GenBank/EMBL | No `molecule_type` annotation (or it was dropped by slicing) | Set `record.annotations['molecule_type'] = 'DNA'` before writing |
| `reverse_complement()`/`transcribe()` returns nonsense, no error | No alphabet validation since 1.78 - a protein/RNA was passed | Track molecule type yourself; only call strand ops on DNA/RNA |

## Decision Tree

```
Need to work with sequence data?
├── Only string-like reads (slice, count, find, translate)?
│   └── Use Seq (immutable)
├── Editing individual positions in place?
│   └── Use MutableSeq, then cast back to Seq to write
├── Need metadata (id, description, features, annotations)?
│   └── Use SeqRecord
└── Writing to GenBank/EMBL?
    └── Use SeqRecord with annotations['molecule_type'] set
```

## Related Skills

- sequence-io/read-sequences - Parse files to get SeqRecord objects
- sequence-io/write-sequences - Write SeqRecord objects to files
- transcription-translation - Transform Seq objects (DNA to protein); inherits the no-alphabet-validation trap
- reverse-complement - Get reverse complement of Seq; silent garbage on non-DNA input
- sequence-slicing - Slice and extract from Seq/SeqRecord; 0-based vs 1-based coordinate trap
- database-access/entrez-fetch - Fetch sequences from NCBI as SeqRecords

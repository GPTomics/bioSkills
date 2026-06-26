---
name: bio-write-sequences
description: Write biological sequences to files (FASTA, FASTQ, GenBank, EMBL) using Biopython Bio.SeqIO. Use when saving sequences, creating new sequence files, or outputting modified records.
tool_type: python
primary_tool: Bio.SeqIO
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Write Sequences

**"Write sequences to a file"** -> Serialize SeqRecord objects into a formatted sequence file.
- Python: `SeqIO.write()` (BioPython)
- R: `writeXStringSet()` (Biostrings)

## Governing Principle

The write is only as complete as the SeqRecord. Each format reads specific record fields and silently ignores the rest, so what survives a write is decided by which fields are populated before the call, not by the format string. FASTA serializes only `id`/`description`+`seq`; FASTQ additionally requires `letter_annotations['phred_quality']`; GenBank/EMBL additionally require `annotations['molecule_type']`. Populate the fields a format needs, or the write either drops data quietly (FASTA) or raises (FASTQ/GenBank).

## Required Import

```python
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
```

## Core Functions

### SeqIO.write() - Write Records to File

```python
SeqIO.write(records, 'output.fasta', 'fasta')
```

- `records` - Single SeqRecord, list, or iterator of SeqRecords
- `handle` - Filename (string) or open file handle
- `format` - Lowercase output format string
- Returns the number of records written (integer)

### record.format() - Get Formatted String

```python
formatted = record.format('fasta')
```

## The FASTA Header Trap (id vs description)

FASTA output is built from `record.description`, NOT `record.id`. The writer compares the first whitespace token of `description` to `id`: if they match it writes `description` as-is; otherwise it prepends `id` + a space. When a record is parsed from FASTA, `description` already leads with the id token, so a round trip is faithful. But when `id` and `description` are set independently, a stale id-like token inside the description gets duplicated, and older BioPython releases dropped the id entirely instead of prepending it.

| record.id | record.description | Header written |
|-----------|--------------------|----------------|
| `seq1` | `seq1 kinase domain` | `>seq1 kinase domain` (clean: description leads with id) |
| `seq1` | `kinase domain` | `>seq1 kinase domain` (id auto-prepended) |
| `seq1` | `` (empty) | `>seq1` (id used as fallback) |
| `seq1` | `gene7 kinase domain` | `>seq1 gene7 kinase domain` (stale id duplicated) |

To control the header exactly and stay robust across versions, make `description` begin with `id` + a space: `description=f'{rec_id} kinase domain'`. The FASTA writer wraps the sequence at 60 characters per line by default.

## Format Field Requirements

| Format | String | Record fields read | Hard requirement |
|--------|--------|--------------------|------------------|
| FASTA | `'fasta'` | id/description, seq | none (header trap above) |
| FASTQ | `'fastq'` | seq, letter_annotations | phred quality scores |
| GenBank | `'genbank'` / `'gb'` | seq, annotations, features | molecule_type |
| EMBL | `'embl'` | seq, annotations, features | molecule_type |
| Tab | `'tab'` | id, seq | none |

## Creating SeqRecord Objects

**Goal:** Construct in-memory records that carry the fields the target format requires.

**Approach:** Build a `SeqRecord` from a `Seq` plus `id`; add `letter_annotations['phred_quality']` for FASTQ and `annotations['molecule_type']` for GenBank/EMBL.

**"Create a sequence record from scratch"** -> Wrap a `Seq` in a `SeqRecord` with metadata.
- Python: `SeqRecord(Seq(...), id=...)` (BioPython)

```python
record = SeqRecord(Seq('ATGCGATCGATCG'), id='seq1', description='seq1 example sequence')
```

## Code Patterns

### Write Single or Multiple Records

```python
records = [SeqRecord(Seq('ATGC'), id='seq1'), SeqRecord(Seq('GCTA'), id='seq2')]
count = SeqIO.write(records, 'output.fasta', 'fasta')
```

### Write to a File Handle (and Append)

```python
with open('output.fasta', 'w') as handle:
    SeqIO.write(records, handle, 'fasta')

with open('output.fasta', 'a') as handle:
    SeqIO.write(new_records, handle, 'fasta')
```

### Write Modified Records via Generator

**Goal:** Transform sequences in memory and write the modified versions to a new file.

**Approach:** Parse input, map a transform over a generator, write the generator. Streaming avoids loading every record into RAM.

**"Modify sequences and save"** -> Parse records, transform each, write with `SeqIO.write()`.

```python
def uppercase_record(rec):
    return SeqRecord(rec.seq.upper(), id=rec.id, description=rec.description)

records = SeqIO.parse('input.fasta', 'fasta')
modified = (uppercase_record(rec) for rec in records)
SeqIO.write(modified, 'output.fasta', 'fasta')
```

### Write FASTQ with Quality Scores

FASTQ requires `letter_annotations['phred_quality']` as a list of ints. `letter_annotations` is length-locked to `len(seq)`: assigning a list whose length differs from the sequence raises. Set the sequence first, then the quality list of matching length.

```python
record = SeqRecord(Seq('ATGCGATCG'), id='read1')
record.letter_annotations['phred_quality'] = [40] * len(record.seq)
SeqIO.write(record, 'output.fastq', 'fastq')
```

### Quality Encoding on Write (Phred vs Solexa)

When both `phred_quality` and `solexa_quality` keys are present, the writer uses Phred. Writing `'fastq-solexa'` from a Phred-only record forces an on-the-fly lossy conversion (the scales diverge in the low-quality region) and emits a `BiopythonWarning` once any score reaches the high end (max quality >= ~62). For modern data, write plain `'fastq'` (Sanger/Phred+33); only use `'fastq-solexa'`/`'fastq-illumina'` when a tool explicitly demands that legacy encoding.

### Write GenBank Format

GenBank and EMBL writing requires `annotations['molecule_type']` (the alphabet that once carried this was removed in BioPython 1.78). Missing it raises on write.

```python
record = SeqRecord(Seq('ATGCGATCGATCG'), id='SEQ001', name='example')
record.annotations['molecule_type'] = 'DNA'
record.annotations['topology'] = 'linear'
record.annotations['organism'] = 'Example organism'
SeqIO.write(record, 'output.gb', 'genbank')
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Header has a duplicated or mangled id | FASTA builds the header from `description`; it does not lead with `id` + space | Set `description=f'{rec.id} ...'` or leave description empty to fall back to id |
| `ValueError: No suitable quality scores found in letter_annotations of SeqRecord (id=...)` on FASTQ write | Record has no `letter_annotations['phred_quality']` | Assign `record.letter_annotations['phred_quality'] = [q]*len(seq)` |
| `TypeError: Any per-letter annotation should be a Python sequence ... of the same length` | Quality list length != `len(seq)` (annotations are length-locked) | Set seq first, then a quality list of matching length |
| `ValueError: missing molecule_type ...` on GenBank/EMBL write | No `annotations['molecule_type']` since the 1.78 alphabet removal | Add `record.annotations['molecule_type'] = 'DNA'` (or 'RNA'/'protein') |
| `BiopythonWarning: Data loss - max Solexa quality ...` | Writing `'fastq-solexa'` from a high Phred-only record forces lossy conversion | Write plain `'fastq'` unless a tool requires the Solexa encoding |
| `TypeError` passing a raw `str`/`Seq` to write | `SeqIO.write` expects SeqRecord(s) | Wrap the sequence in a `SeqRecord` first |
| `ValueError: Sequences must all be the same length` | PHYLIP/alignment format with unequal lengths | Align, pad, or trim to equal length first |

## Related Skills

- read-sequences - Read sequences before modifying and writing
- format-conversion - Direct format conversion without intermediate processing
- filter-sequences - Filter sequences before writing a subset
- fastq-quality - Phred/Solexa encodings and quality-score handling
- sequence-manipulation/seq-objects - Create SeqRecord objects to write
- alignment-files/sam-bam-basics - For SAM/BAM output, use samtools/pysam

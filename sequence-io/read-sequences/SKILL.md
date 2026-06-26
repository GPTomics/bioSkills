---
name: bio-read-sequences
description: Read biological sequence files (FASTA, FASTQ, GenBank, EMBL, ABI, SFF) with Biopython Bio.SeqIO, choosing between streaming, in-memory, and on-disk-indexed access. Use when parsing sequence files, iterating multi-record files, randomly accessing records by ID in large files, or maximizing parse throughput.
tool_type: python
primary_tool: Bio.SeqIO
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show biopython` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Read Sequences

Read biological sequence data from files using Biopython's Bio.SeqIO module.

**"Read sequences from a file"** -> Parse a file into SeqRecord objects exposing id, sequence, and annotations.
- Python: `SeqIO.parse()` / `SeqIO.read()` (BioPython)
- R: `readDNAStringSet()` / `readAAStringSet()` (Biostrings)

## The Governing Principle

Stream by default. `SeqIO.parse()` yields one record at a time and never holds the whole file in RAM, so it scales to any size. Reach for an in-memory or indexed structure only when the access pattern demands it: load all records (`to_dict`) only for small files needing random access; build an index (`index` / `index_db`) for random access into large files. Never `list()` a huge file or `to_dict()` it - that defeats streaming and can exhaust memory.

## Which Function to Use

| Method | Returns | Memory model | Random access | Persists | Multi-file |
|--------|---------|--------------|---------------|----------|------------|
| `parse(handle, format)` | generator of SeqRecord | one record at a time | no | no | no |
| `read(handle, format)` | one SeqRecord | one record | n/a | no | no |
| `to_dict(records)` | real `dict` | ALL records in RAM | yes | no | feed combined iterators |
| `index(filename, format)` | dict-like (read-only) | byte offsets only, re-parses on access | yes | no | no |
| `index_db(idx_file, files, format)` | dict-like (read-only) | on-disk SQLite index | yes | yes | yes |

Decision rule: `parse` for streaming; `read` for a known single-record file; `to_dict` when the file is small and random access by ID is needed; `index` for random access into one large file; `index_db` for files larger than RAM, many files indexed together, or an index reused across runs.

Behavioral traps these methods hide:
- `parse()` is a one-pass generator. It is NOT subscriptable (`parse(...)[3]` raises TypeError), and it EXHAUSTS SILENTLY: a second `for` loop over the same generator object yields nothing with no error. Re-call `parse()` for each pass, or `list()` it once if the file is small.
- `read()` fails LOUDLY: zero records raise `ValueError: No records found in handle`; more than one raises `ValueError: More than one record found in handle`. Use it as an assertion that the file holds exactly one sequence.
- `to_dict()`, `index()`, and `index_db()` all raise `ValueError` on a DUPLICATE id (`Duplicate key '...'`). Supply a `key_function` to derive unique keys when ids collide.
- `index()` needs a FILENAME, not a handle (it must seek). It stores only byte offsets and re-parses the record from disk on every access, so it returns a fresh object each time and mutations do not persist. It is read-only (`__setitem__` raises NotImplementedError).
- `index_db()` stores the offset index in an on-disk SQLite file. It PERSISTS across sessions (reopen later with just the index filename), and scales beyond RAM and across multiple files (pass a list of filenames). This is the right answer for data larger than memory.

The `alphabet=` argument still appears in some signatures for back-compatibility but is a no-op since BioPython 1.78; leave it `None`.

## Required Import

```python
from Bio import SeqIO
```

## Reading Records

### SeqIO.parse() - Stream Multiple Records
Returns a one-pass iterator of SeqRecord objects. Always pass the format explicitly as the second argument.

```python
for record in SeqIO.parse('sequences.fasta', 'fasta'):
    print(record.id, len(record.seq))
```

### SeqIO.read() - Exactly One Record
Use when the file must contain a single sequence; raises on zero or multiple records.

```python
record = SeqIO.read('single.fasta', 'fasta')
```

## Random Access

### SeqIO.to_dict() - Small Files
Loads every record into a dictionary keyed by id. Fast random access, but holds all records in RAM.

```python
records = SeqIO.to_dict(SeqIO.parse('sequences.fasta', 'fasta'))
seq = records['sequence_id'].seq
```

### SeqIO.index() - One Large File

**Goal:** Random access by id into a large file without loading every record into memory.

**Approach:** Build an in-memory map of byte offsets keyed by id; each lookup re-parses one record from disk.

**Reference (BioPython 1.83+):**
```python
records = SeqIO.index('large.fasta', 'fasta')
seq = records['sequence_id'].seq
records.close()
```

A `key_function` maps the id STRING to a custom key (note: `to_dict`'s key_function receives the whole record instead):
```python
def get_accession(identifier):
    return identifier.split('.')[0]  # drop the version suffix

records = SeqIO.index('sequences.fasta', 'fasta', key_function=get_accession)
```

### SeqIO.index_db() - Huge / Multiple Files

**Goal:** Random access into data larger than RAM, or across many files, with the index reusable across runs.

**Approach:** Persist the offset index in an on-disk SQLite database; reopen it later without re-parsing.

**Reference (BioPython 1.83+):**
```python
# First call parses the file(s) and builds the SQLite index
records = SeqIO.index_db('index.sqlite', 'large.fasta', 'fasta')
seq = records['sequence_id'].seq
records.close()

# Later sessions reopen instantly with just the index filename
records = SeqIO.index_db('index.sqlite')

# Index multiple files as one database
records = SeqIO.index_db('combined.sqlite', ['file1.fasta', 'file2.fasta'], 'fasta')
```

## High-Performance Parsing

For maximum throughput on large files, low-level parsers (SimpleFastaParser, FastqGeneralIterator) yield raw tuples and skip SeqRecord construction, so they run substantially faster than SeqIO.parse.

### SimpleFastaParser

**Goal:** Parse large FASTA files at maximum speed without SeqRecord overhead.

**Approach:** Iterate `(title, sequence)` string tuples directly from the handle.

**Reference (BioPython 1.83+):**
```python
from Bio.SeqIO.FastaIO import SimpleFastaParser

with open('large.fasta') as handle:
    for title, sequence in SimpleFastaParser(handle):
        if len(sequence) > 1000:
            seq_id = title.split()[0]  # first whitespace token is the id
```

### FastqGeneralIterator

**Goal:** Parse large FASTQ files at maximum speed.

**Approach:** Iterate `(title, sequence, quality_string)` string tuples; decode quality manually if needed.

**Reference (BioPython 1.83+):**
```python
from Bio.SeqIO.QualityIO import FastqGeneralIterator

with open('reads.fastq') as handle:
    for title, sequence, quality in FastqGeneralIterator(handle):
        avg_qual = sum(ord(c) - 33 for c in quality) / len(quality)  # Phred+33
```

## SeqRecord Attributes

After parsing, each record exposes:

```python
record.id          # first whitespace token of the header (string)
record.name        # same first token (for FASTA, name == id)
record.description # the ENTIRE header after '>', including the id token
record.seq         # sequence data (Seq object; case-preserving)
record.features    # list of SeqFeature objects (GenBank/EMBL)
record.annotations # dict of annotations (organism, molecule_type, ...)
record.letter_annotations  # per-letter dict (e.g. 'phred_quality' list)
record.dbxrefs     # database cross-references
```

### id vs name vs description - the first-space split
A FASTA header `>FIRST rest of the line` parses to: `id` = `FIRST` (the first whitespace token), `name` = `FIRST` (same token), `description` = `FIRST rest of the line` (the WHOLE header after `>`, including the id). So `>seq1 some desc` gives id `seq1`, name `seq1`, description `seq1 some desc`. The id is therefore the leading word of the description, not a separate field - relevant when writing records back out.

## Common Formats

| Format | String | Typical Extension | Notes |
|--------|--------|-------------------|-------|
| FASTA | `'fasta'` | .fasta, .fa, .fna, .faa | Most common |
| FASTA 2-line | `'fasta-2line'` | .fasta | One line per sequence (no wrapping) |
| FASTQ | `'fastq'` | .fastq, .fq | Alias of fastq-sanger (Phred+33) |
| FASTQ Solexa | `'fastq-solexa'` | .fastq | Old Solexa (Solexa+64, scores -5..62) |
| FASTQ Illumina | `'fastq-illumina'` | .fastq | Illumina 1.3-1.7 (Phred+64) |
| GenBank | `'genbank'` or `'gb'` | .gb, .gbk | With features/annotations |
| EMBL | `'embl'` | .embl | European format with features |
| Swiss-Prot | `'swiss'` | .dat | UniProt format |

FASTQ quality encoding cannot be auto-detected reliably: the same quality line can be valid Phred+33 and Phred+64. Picking the wrong string can silently shift every score by 31. Confirm the encoding before parsing; see fastq-quality for the full encoding decision.

## Specialized Formats

| Format | String | Use Case |
|--------|--------|----------|
| ABI | `'abi'` | Sanger sequencing trace files (.ab1) |
| ABI Trimmed | `'abi-trim'` | ABI with low-quality ends trimmed |
| SFF | `'sff'` | 454/Ion Torrent flowgram data |
| SFF Trimmed | `'sff-trim'` | SFF with adapter/quality trimming |
| QUAL | `'qual'` | Quality scores file (pairs with FASTA) |
| PDB SEQRES | `'pdb-seqres'` | Protein sequences from PDB SEQRES records |
| PDB ATOM | `'pdb-atom'` | Sequences from ATOM records in PDB |
| SnapGene | `'snapgene'` | SnapGene .dna files |

### Reading ABI Trace Files
```python
record = SeqIO.read('sample.ab1', 'abi')
qualities = record.letter_annotations['phred_quality']
record_trimmed = SeqIO.read('sample.ab1', 'abi-trim')  # low-quality ends removed
```

### Reading 454/Ion Torrent SFF
```python
for record in SeqIO.parse('reads.sff', 'sff'):
    print(record.id, len(record.seq))
```

### Reading PDB Sequences
```python
for record in SeqIO.parse('structure.pdb', 'pdb-seqres'):
    print(record.id, record.seq)
```

## Alignment Formats (Read-Only)

| Format | String | Notes |
|--------|--------|-------|
| PHYLIP | `'phylip'` | Interleaved; `'phylip-relaxed'` allows longer names |
| Clustal | `'clustal'` | ClustalW output |
| Stockholm | `'stockholm'` | Rfam/Pfam alignments |
| NEXUS | `'nexus'` | PAUP/MrBayes format |
| MAF | `'maf'` | Multiple Alignment Format |

## Code Patterns

### Count Records Without Loading All
```python
count = sum(1 for _ in SeqIO.parse('sequences.fasta', 'fasta'))
```

### Read GenBank with Features
```python
for record in SeqIO.parse('sequence.gb', 'genbank'):
    for feature in record.features:
        if feature.type == 'CDS':
            product = feature.qualifiers.get('product', ['Unknown'])[0]
            cds_seq = feature.extract(record.seq)  # spliced feature sequence
```

### Access FASTQ Quality Scores
```python
for record in SeqIO.parse('reads.fastq', 'fastq'):
    qualities = record.letter_annotations['phred_quality']
    avg_quality = sum(qualities) / len(qualities)
```

### Read From a File Handle
```python
with open('sequences.fasta') as handle:
    for record in SeqIO.parse(handle, 'fasta'):
        print(record.id)
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Second loop over a parser yields nothing, no error | `parse()` generator exhausted after the first pass | Re-call `parse()` per pass, or `list()` once for small files |
| `TypeError: 'generator' object is not subscriptable` | Indexed/sliced a `parse()` result | Wrap in `list()`, or use `to_dict`/`index` for keyed access |
| `ValueError: More than one record found in handle` | `read()` on a multi-record file | Use `parse()` |
| `ValueError: No records found in handle` | `read()` on an empty/zero-record file | Check the file and format string; use `parse()` if multi-record |
| `ValueError: Duplicate key '...'` | `to_dict`/`index`/`index_db` hit a repeated id | Pass a `key_function` that derives unique keys |
| Random access by id silently slow / re-reads disk | `index()` re-parses each access; mutations don't persist | Expected; cache needed records, or use `to_dict` for small files |
| MemoryError / process killed on a huge file | `list()` or `to_dict()` loaded everything into RAM | Stream with `parse()`; use `index_db()` for random access |
| `ValueError: unknown format` | Misspelled format string | Use a lowercase string from the format tables |
| `ValueError`/`AssertionError` naming the LOCUS line | GenBank parser reads fixed LOCUS columns (molecule type ~44-54, topology ~55-63); ICE/SnapGene/Ensembl/assembler LOCUS lines violate the spec | Biologically valid content can still fail the strict column parse; fix the LOCUS columns or re-export from a spec-compliant writer |
| FASTQ scores all off by ~31 with no error | Wrong FASTQ variant string (Phred+33 vs +64 overlap) | Confirm encoding; see fastq-quality |
| `AttributeError` referencing `.alphabet` | Code assumes pre-1.78 alphabet API | Drop alphabet usage; molecule type lives in `annotations['molecule_type']` |

## Related Skills

- write-sequences - Write parsed sequences to new files
- filter-sequences - Filter sequences by criteria after reading
- format-conversion - Convert between formats (GenBank->FASTA silently drops annotations)
- compressed-files - Read gzip/bzip2/BGZF compressed files; only BGZF supports indexed random access
- fastq-quality - FASTQ encoding (Phred vs Solexa) and offset selection
- sequence-manipulation/seq-objects - Work with parsed SeqRecord and Seq objects
- database-access/entrez-fetch - Fetch sequences from NCBI instead of local files
- alignment-files/sam-bam-basics - For SAM/BAM/CRAM alignment files, use samtools/pysam

---
name: bio-filter-sequences
description: Filter and select sequences by criteria (length, ID, GC content, N content, motifs, patterns, description) using Biopython, streaming so large files never load into RAM. Use when subsetting a FASTA/FASTQ file, removing unwanted or low-quality records, or selecting records by specific criteria. Use the paired-end-fastq skill instead whenever the input is paired R1/R2 reads.
tool_type: python
primary_tool: Bio.SeqIO
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Filter Sequences

**"Filter sequences by length, quality, or content"** -> Apply boolean criteria to a stream of sequence records and write survivors to output.
- Python: generator expression with `SeqIO.parse()` + `SeqIO.write()` (BioPython)
- CLI: `seqkit seq -m 200` (SeqKit) or `awk` on FASTA

Filter and select sequences based on length, ID, GC content, N content, motifs, regex patterns, and description.

## The governing principle

Two traps cause silent, downstream-corrupting errors. Neither raises an exception, so the agent must guard against both up front.

1. **Never filter one mate of a paired-end set independently.** Aligners (bwa, bowtie2) read R1 and R2 as two parallel streams and pair the i-th record of each, assuming SAME ORDER and SAME COUNT. Dropping a read from R1 alone DESYNCS the files: best case the aligner crashes on a name mismatch, worst case it silently pairs the wrong R1 with the wrong R2, producing mismapped reads and garbage insert sizes with no error. If the input is paired, route to the paired-end-fastq skill (synchronized filtering that writes matched output plus separate orphan/singleton files). Do NOT apply the single-file patterns below to one mate.
2. **Filter by STREAMING, not by loading.** `SeqIO.parse()` yields one record at a time; a generator expression into `SeqIO.write()` holds a single record in RAM regardless of file size. `list(SeqIO.parse(...))` materializes every record and OOMs on large FASTQ. Only the patterns that genuinely need all records at once (random sampling, splitting into multiple files) load the file; they say so explicitly.

## Required Imports

```python
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
```

## Core Pattern

Stream records through a generator expression so memory stays flat:

```python
records = SeqIO.parse('input.fasta', 'fasta')
filtered = (rec for rec in records if len(rec.seq) >= 100)
SeqIO.write(filtered, 'output.fasta', 'fasta')
```

`SeqIO.write()` consumes the generator lazily and returns the count written.

## Filter by Length

### Minimum Length
```python
records = SeqIO.parse('input.fasta', 'fasta')
long_seqs = (rec for rec in records if len(rec.seq) >= 500)
SeqIO.write(long_seqs, 'long.fasta', 'fasta')
```

### Length Range
```python
records = SeqIO.parse('input.fasta', 'fasta')
sized = (rec for rec in records if 100 <= len(rec.seq) <= 1000)
SeqIO.write(sized, 'sized.fasta', 'fasta')
```

### Remove Short Sequences
```python
min_length = 200
records = SeqIO.parse('input.fasta', 'fasta')
filtered = (rec for rec in records if len(rec.seq) >= min_length)
count = SeqIO.write(filtered, 'filtered.fasta', 'fasta')
```

`len(rec.seq)` counts every base, including soft-masked lowercase (see Case Sensitivity below).

## Filter by ID

### Select Specific IDs
```python
wanted_ids = {'seq1', 'seq2', 'seq3'}
records = SeqIO.parse('input.fasta', 'fasta')
selected = (rec for rec in records if rec.id in wanted_ids)
SeqIO.write(selected, 'selected.fasta', 'fasta')
```

`rec.id` is the first whitespace-delimited token of the header. For CASAVA 1.8+ paired FASTQ, R1 and R2 share this id (the mate number lives after the space), so an id set matches both mates equally - another reason mate-aware filtering belongs in paired-end-fastq.

### Select from ID File

**Goal:** Extract sequences whose IDs appear in an external list file.

**Approach:** Load IDs into a set for O(1) lookup, then stream-filter and write matches.

**Reference (BioPython 1.83+):**
```python
with open('ids.txt') as f:
    wanted_ids = {line.strip() for line in f}

records = SeqIO.parse('input.fasta', 'fasta')
selected = (rec for rec in records if rec.id in wanted_ids)
SeqIO.write(selected, 'selected.fasta', 'fasta')
```

### Exclude Specific IDs
```python
exclude_ids = {'bad_seq1', 'bad_seq2'}
records = SeqIO.parse('input.fasta', 'fasta')
kept = (rec for rec in records if rec.id not in exclude_ids)
SeqIO.write(kept, 'kept.fasta', 'fasta')
```

### Filter by ID Pattern
```python
import re

pattern = re.compile(r'^chr\d+$')  # matches chr1, chr2, etc.
records = SeqIO.parse('input.fasta', 'fasta')
chromosomes = (rec for rec in records if pattern.match(rec.id))
SeqIO.write(chromosomes, 'chromosomes.fasta', 'fasta')
```

## Filter by GC Content

**Goal:** Keep records whose GC fraction falls in a target band.

**Approach:** Use `gc_fraction()`, which returns a FRACTION (0-1), NOT a percentage - thresholds must be 0.4, not 40. The `ambiguous=` mode decides how N and other IUPAC ambiguity codes are counted, and the same sequence yields a different GC value per mode, so set it explicitly rather than relying on the default.

**Reference (BioPython 1.83+):**
```python
from Bio.SeqUtils import gc_fraction

records = SeqIO.parse('input.fasta', 'fasta')
moderate_gc = (rec for rec in records if 0.4 <= gc_fraction(rec.seq, ambiguous='ignore') <= 0.6)
SeqIO.write(moderate_gc, 'moderate_gc.fasta', 'fasta')
```

### Choosing the ambiguous= mode

`gc_fraction(seq, ambiguous='remove')` is the default. For the same sequence the three modes give different answers - an N-containing read can pass or fail purely because of the mode:

| Mode | Denominator | `gc_fraction('GCGCNNNN')` | When to use |
|------|-------------|---------------------------|-------------|
| `'remove'` (default) | only unambiguous A,T,G,C,S,W,U | 1.0 | GC of the called bases only; ignores how many N's are present |
| `'ignore'` | full `len(seq)` (N's dilute GC) | 0.5 | GC over the whole read; matches a naive `(G+C)/len` |
| `'weighted'` | full length, ambiguous codes add expected GC | 0.75 | each IUPAC code contributes its mean GC (S=1.0, W=0.0, N=0.5, V/B=0.667, H/D=0.333) |

A naive `(G+C)/len` silently equals `'ignore'` mode and under-reports GC whenever N's are present. The default `'remove'` ignores N's entirely, so a heavily-N read can post a misleadingly extreme GC. Pick the mode that matches the intent and pass it explicitly.

### High / Low GC bands
```python
records = SeqIO.parse('input.fasta', 'fasta')
high_gc = (rec for rec in records if gc_fraction(rec.seq, ambiguous='ignore') >= 0.6)
SeqIO.write(high_gc, 'high_gc.fasta', 'fasta')
```

## Case Sensitivity (soft-masking)

`Seq` is CASE-PRESERVING: lowercase soft-masked bases (from RepeatMasker, Ensembl, dustmasker) survive parse and round-trip unchanged. Length, GC, motif, and regex filters are CASE-SENSITIVE - a naive uppercase test silently misses masked bases. Always `.upper()` the sequence before content matching when the masking should not affect the decision:

```python
seq_upper = str(rec.seq).upper()
has_site = 'GAATTC' in seq_upper            # matches gaattc and GAATTC
```

`gc_fraction()` itself is case-insensitive, but a hand-rolled `.count('G')` is not - count on the uppercased string.

## Filter by Sequence Content

### Remove Sequences with N's
```python
records = SeqIO.parse('input.fasta', 'fasta')
clean = (rec for rec in records if 'N' not in str(rec.seq).upper())
SeqIO.write(clean, 'clean.fasta', 'fasta')
```

### Limit N Content
```python
def n_fraction(seq):
    upper = str(seq).upper()
    return upper.count('N') / len(seq)

records = SeqIO.parse('input.fasta', 'fasta')
low_n = (rec for rec in records if n_fraction(rec.seq) < 0.05)  # under 5% ambiguous bases
SeqIO.write(low_n, 'low_n.fasta', 'fasta')
```

### Contains Specific Motif
```python
motif = 'GAATTC'  # EcoRI site
records = SeqIO.parse('input.fasta', 'fasta')
with_motif = (rec for rec in records if motif in str(rec.seq).upper())
SeqIO.write(with_motif, 'with_ecori.fasta', 'fasta')
```

### Regex Pattern in Sequence
```python
import re

pattern = re.compile(r'ATG.{30,100}T(AA|AG|GA)')  # ORF-like pattern
records = SeqIO.parse('input.fasta', 'fasta')
matches = (rec for rec in records if pattern.search(str(rec.seq).upper()))
SeqIO.write(matches, 'orf_like.fasta', 'fasta')
```

## Filter by Description

### Description Contains Keyword
```python
records = SeqIO.parse('input.fasta', 'fasta')
kinases = (rec for rec in records if 'kinase' in rec.description.lower())
SeqIO.write(kinases, 'kinases.fasta', 'fasta')
```

### Multiple Keywords (OR)
```python
keywords = ['kinase', 'phosphatase', 'transferase']
records = SeqIO.parse('input.fasta', 'fasta')
enzymes = (rec for rec in records if any(k in rec.description.lower() for k in keywords))
SeqIO.write(enzymes, 'enzymes.fasta', 'fasta')
```

## Combine Multiple Filters

**Goal:** Remove sequences that fail any of several length/content thresholds.

**Approach:** Define a predicate that checks all criteria against the uppercased sequence once, set the GC `ambiguous=` mode explicitly, apply the predicate as a generator filter, and stream survivors to output.

**Reference (BioPython 1.83+):**
```python
from Bio.SeqUtils import gc_fraction

def passes_filters(record):
    if len(record.seq) < 100:
        return False
    gc = gc_fraction(record.seq, ambiguous='ignore')
    if gc < 0.3 or gc > 0.7:
        return False
    if 'N' in str(record.seq).upper():
        return False
    return True

records = SeqIO.parse('input.fasta', 'fasta')
filtered = (rec for rec in records if passes_filters(rec))
SeqIO.write(filtered, 'filtered.fasta', 'fasta')
```

## Sample Sequences

### Random Sample (requires loading all)
```python
import random

records = list(SeqIO.parse('input.fasta', 'fasta'))  # loads file - needs all records up front
sample = random.sample(records, min(100, len(records)))
SeqIO.write(sample, 'sample.fasta', 'fasta')
```

### First N Sequences (streams)
```python
from itertools import islice

records = SeqIO.parse('input.fasta', 'fasta')
first_100 = islice(records, 100)
SeqIO.write(first_100, 'first100.fasta', 'fasta')
```

### Every Nth Sequence (streams)
```python
records = SeqIO.parse('input.fasta', 'fasta')
every_10th = (rec for i, rec in enumerate(records) if i % 10 == 0)
SeqIO.write(every_10th, 'sampled.fasta', 'fasta')
```

## Split by Criteria

### Split by Length

**Goal:** Partition sequences into separate files based on a length threshold.

**Approach:** Load all records once, partition with list comprehensions, and write each partition. Loading is acceptable here because both partitions are needed in a single pass; for very large files, run two streaming passes instead.

**Reference (BioPython 1.83+):**
```python
records = list(SeqIO.parse('input.fasta', 'fasta'))
short = [r for r in records if len(r.seq) < 500]
long = [r for r in records if len(r.seq) >= 500]
SeqIO.write(short, 'short.fasta', 'fasta')
SeqIO.write(long, 'long.fasta', 'fasta')
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Downstream mismapping, wrong insert sizes, no error | Filtered one mate of a paired-end set independently, desyncing R1/R2 | Never filter one mate alone; use paired-end-fastq for synchronized filtering with orphan output |
| GC filter keeps/drops the wrong reads | `gc_fraction` returns a fraction 0-1 but threshold written as a percent (40 instead of 0.4) | Use 0-1 thresholds; multiply by 100 only for display |
| N-containing read unexpectedly passes or fails GC band | Wrong `ambiguous=` mode (default `'remove'` drops N's; `'ignore'` dilutes GC) | Set `ambiguous=` explicitly to match intent |
| Soft-masked read fails a motif/regex/uppercase test | `Seq` is case-preserving; lowercase masked bases do not match an uppercase pattern | `.upper()` the sequence before content matching |
| Generator yields nothing on second use | `SeqIO.parse()` is one-pass and exhausts silently | Re-create the generator, or `list()` it if it must be reused |
| MemoryError on large FASTQ | `list(SeqIO.parse(...))` materialized every record | Use a generator expression; only load for sampling/splitting |
| Empty output file | Filter too strict, or matched against the wrong case/field | Loosen thresholds; confirm id vs description and case |

## Related Skills

- read-sequences - Parse sequences before filtering
- write-sequences - Write filtered sequences to output
- fastq-quality - Filter FASTQ by per-base quality scores and encoding
- paired-end-fastq - Synchronized filtering of R1/R2 with orphan handling
- sequence-manipulation/sequence-properties - Per-sequence GC, length, and composition
- sequence-manipulation/motif-search - Filter by complex motif patterns
- alignment-files/alignment-filtering - Filter aligned reads with samtools view -f/-F

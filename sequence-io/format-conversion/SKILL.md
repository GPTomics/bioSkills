---
name: bio-format-conversion
description: Convert between sequence file formats (FASTA, FASTQ, GenBank, EMBL, Stockholm) and re-encode FASTQ quality offsets using Biopython Bio.SeqIO. Use when changing a file format for a downstream tool, fixing FASTQ quality encoding (Phred+33 vs Phred+64 vs Solexa), or when a conversion risks silently dropping annotations or quality scores.
tool_type: python
primary_tool: Bio.SeqIO
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Format Conversion

**"Convert this file to a different format"** -> Read records in one format, optionally add or drop annotations, and write in the target format.
- Python: `SeqIO.convert()` for a streaming one-shot conversion, or `SeqIO.parse()` + `SeqIO.write()` when records need modification (BioPython)
- CLI: `seqkit seq` (SeqKit) for FASTA/FASTQ; `samtools view` for SAM/BAM/CRAM

## The Governing Principle

A conversion is lossy whenever the target format cannot represent the source's information. The conversion still succeeds with no error and no warning. FASTA stores only id + description + sequence, so it is the most lossy common target: converting GenBank, EMBL, or FASTQ to FASTA silently discards everything the richer format carried. Before converting, decide whether the destination can hold what the source contains. If it cannot, treat the conversion as a deliberate downgrade, not a neutral reformat.

## The Canonical Trap: GenBank/EMBL -> FASTA Silently Drops Everything

`SeqIO.convert('in.gb', 'genbank', 'out.fasta', 'fasta')` discards all features, annotations, qualifiers, and dbxrefs. The genes, CDS coordinates, `/product` and `/gene` qualifiers, organism, taxonomy, references, and molecule_type are all gone. There is no error, no warning, and the record count is unchanged, so the loss is invisible unless the output is inspected. FASTA encodes only `record.id`, `record.description`, and `record.seq`; everything in `record.features`, `record.annotations`, and `record.dbxrefs` has nowhere to go.

If the features matter, do not convert to FASTA. Extract feature sequences first (see sequence-manipulation/sequence-slicing) or keep the GenBank file as the source of record and use the FASTA only as a sequence-only derivative for tools that demand FASTA.

## Lossy Conversion Decision Table

| From | To | What is lost (silently) |
|------|-----|-------------------------|
| GenBank / EMBL | FASTA | All features, qualifiers, annotations, dbxrefs; keeps id + description + seq |
| GenBank | EMBL (or reverse) | Usually lossless; both hold features and annotations |
| FASTQ | FASTA | Per-base quality scores (`phred_quality`) |
| FASTQ Phred+64 | FASTQ Phred+33 | Nothing if offsets handled correctly; corruption if the wrong parser is used |
| FASTQ Phred | FASTQ Solexa | Precision at low quality (round-trip lossy below ~Q10); warns when max Solexa exceeded |
| Stockholm | FASTA | Alignment columns (gaps), consensus, per-column annotation; keeps ungapped seqs |
| Any rich format | FASTA | Everything except id + description + seq |

The general pattern: rich -> flat loses the richness. The conversion succeeds regardless.

## Preferred Path: SeqIO.convert() Is a Streaming One-Shot

For a plain conversion with no record modification, use `SeqIO.convert()`. It streams one record at a time from input to output (memory-efficient, never loads the whole file) and is preferred over `parse()` + `write()`, which is only needed when records must be changed en route.

```python
from Bio import SeqIO

count = SeqIO.convert('input.gb', 'genbank', 'output.fasta', 'fasta')
print(f'Converted {count} records')
```

Parameters: `in_file`, `in_format`, `out_file`, `out_format` (filenames or handles; format strings are lowercase). Returns the number of records written. Reach for `parse()` + `write()` only when injecting annotations, transforming sequences, or filtering during the conversion.

## FASTQ Quality Encoding Conversion

FASTQ quality is one ASCII character per base, but the offset and score type differ across instrument generations. Re-encoding between them is a conversion, not a copy: the bytes in the quality line change.

| Format string | For | Offset | Score type |
|---------------|-----|--------|-----------|
| `fastq` (alias of `fastq-sanger`) | Sanger and modern Illumina 1.8+ | 33 | Phred 0-93 |
| `fastq-sanger` | same as above | 33 | Phred 0-93 |
| `fastq-illumina` | Illumina 1.3-1.7 | 64 | Phred 0-62 |
| `fastq-solexa` | pre-1.3 Solexa | 64 | Solexa odds -5..62 |

Re-encode old Illumina 1.3+ (Phred+64) to modern Sanger (Phred+33) by naming both variants. `SeqIO.convert()` reads with the input offset and writes with the output offset:

```python
from Bio import SeqIO

SeqIO.convert('illumina13.fastq', 'fastq-illumina', 'sanger.fastq', 'fastq-sanger')
```

**Never re-encode without a verified source encoding.** Quality encoding cannot be auto-detected in general: ASCII >= 64 is legal in every variant, so a high-quality Sanger file and a low-quality Illumina-1.3 file can be byte-identical in their quality lines. Two failure modes follow from guessing wrong:
- Loud and safe: a character outside the chosen parser's range raises `ValueError` noting the quality string is not in the correct range for the chosen QualityIO parser.
- Silent and dangerous: if every quality char falls in the overlap valid for both encodings, no error fires and every score comes out off by exactly 31 (the 64 - 33 offset gap). QC is then silently garbage.

Confirm the encoding from the sequencing pipeline (or FastQC's inferred encoding) before re-encoding; do not let the agent guess the offset.

**Solexa is doubly lossy.** Solexa uses an odds score, Q = -10 log10(P/(1-P)), not Phred's Q = -10 log10(P), which is why Solexa scores go negative. Phred <-> Solexa conversions round a float to one ASCII char per base, so the round trip is many-to-one and lossy below ~Q10 (for example Solexa 9 and 10 both map to Phred 10). Writing `fastq-solexa` from a Phred-only record forces an on-the-fly lossy conversion and emits a `BiopythonWarning` when `max(qualities) >= 62.5`. There is no clean Phred -> Solexa path that avoids the loss; only re-encode toward Solexa when a legacy tool truly requires it.

## Conversions That Require Adding Data

FASTA has no molecule_type and no quality, so converting FASTA up to a richer format means supplying what FASTA lacked. Stream records through a generator that injects the missing field.

### FASTA to GenBank (requires molecule_type)

**Goal:** Convert FASTA to GenBank, which the writer refuses to produce without `molecule_type`.

**Approach:** Stream records through a generator that sets `record.annotations['molecule_type']`, then write as GenBank.

**Reference (BioPython 1.83+):**
```python
from Bio import SeqIO

def add_molecule_type(records, mol_type='DNA'):
    for record in records:
        record.annotations['molecule_type'] = mol_type
        yield record

records = SeqIO.parse('input.fasta', 'fasta')
SeqIO.write(add_molecule_type(records), 'output.gb', 'genbank')
```

### FASTA to FASTQ (requires quality scores)

**Goal:** Convert FASTA to FASTQ by assigning placeholder per-base quality.

**Approach:** Stream records through a generator that adds a `phred_quality` list of the right length, then write as FASTQ.

**Reference (BioPython 1.83+):**
```python
from Bio import SeqIO

def add_quality(records, quality=40):
    for record in records:
        record.letter_annotations['phred_quality'] = [quality] * len(record.seq)
        yield record

records = SeqIO.parse('input.fasta', 'fasta')
SeqIO.write(add_quality(records), 'output.fastq', 'fastq')
```

Placeholder quality is fabricated, not measured: downstream QC and variant callers will treat it as real. Use it only to satisfy a tool's format requirement, never to imply the bases were measured at that quality. `letter_annotations` is length-locked to the sequence, so the list length must equal `len(record.seq)`.

## Batch Convert a Directory

**Goal:** Convert every file of one format in a directory to another format.

**Approach:** Glob the input files, apply `SeqIO.convert()` to each, and report per-file counts.

**Reference (BioPython 1.83+):**
```python
from pathlib import Path
from Bio import SeqIO

for gb_file in Path('.').glob('*.gb'):
    fasta_file = gb_file.with_suffix('.fasta')
    count = SeqIO.convert(str(gb_file), 'genbank', str(fasta_file), 'fasta')
    print(f'{gb_file.name}: {count} records')
```

## Convert With Sequence Modification

When the conversion must also transform sequences, parse and write explicitly rather than using `convert()`.

```python
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def uppercase_record(rec):
    return SeqRecord(rec.seq.upper(), id=rec.id, description=rec.description)

records = SeqIO.parse('input.fasta', 'fasta')
SeqIO.write((uppercase_record(rec) for rec in records), 'output.fasta', 'fasta')
```

`Seq` is case-preserving, so lowercase soft-masking survives a plain conversion; call `.upper()` explicitly only when the destination tool requires uppercase.

## Alignment Format Conversion

Sequence formats drop gaps and alignment columns. To convert between alignment formats (Stockholm, PHYLIP, Clustal, FASTA-alignment) keeping the columns, use `AlignIO`, not `SeqIO`.

```python
from Bio import AlignIO

AlignIO.convert('alignment.sto', 'stockholm', 'alignment.phy', 'phylip')
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| GenBank features missing after conversion | Target was FASTA, which cannot hold features | Expected and silent; keep the GenBank as source, or extract features before converting |
| `ValueError` about missing molecule_type | Writing GenBank/EMBL from records that lack it (e.g. from FASTA) | Set `record.annotations['molecule_type']` before writing |
| `ValueError` about quality scores | Writing FASTQ from records with no `phred_quality` | Add `phred_quality` to `letter_annotations` (length must equal the sequence) |
| `ValueError` mentioning the QualityIO parser | A quality char is outside the named parser's range (wrong FASTQ variant) | Use the correct variant: `fastq-sanger`, `fastq-illumina`, or `fastq-solexa` |
| FASTQ scores all off by 31 with no error | Read Phred+33 as `fastq-illumina` or Phred+64 as `fastq-sanger` (overlap region) | Confirm the true encoding from the pipeline; re-read with the right variant |
| `BiopythonWarning` "Data loss - max Solexa quality" | Writing `fastq-solexa` from Phred scores above ~62 | Expected lossy conversion; only write Solexa when a legacy tool requires it |
| Alignment columns/gaps lost | Used `SeqIO` on an alignment | Use `AlignIO.convert()` to preserve columns |

## Related Skills

- read-sequences - Parse sequences and choose parse vs index for the input
- write-sequences - Write converted sequences with modifications
- fastq-quality - Phred/Solexa/Illumina encoding details and quality handling
- batch-processing - Convert many files across a directory
- compressed-files - Handle gzip/BGZF input and output during conversion
- sequence-manipulation/sequence-slicing - Extract feature sequences before downgrading to FASTA
- alignment-files/sam-bam-basics - For SAM/BAM/CRAM conversion, use samtools view

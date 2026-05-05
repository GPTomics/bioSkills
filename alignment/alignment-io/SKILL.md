---
name: bio-alignment-io
description: Read, write, and convert multiple sequence alignment files using Biopython Bio.AlignIO. Supports Clustal, PHYLIP, Stockholm, FASTA, Nexus, and other alignment formats for phylogenetics and conservation analysis. Use when reading, writing, or converting alignment file formats.
tool_type: python
primary_tool: Bio.AlignIO
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Alignment File I/O

Read, write, and convert multiple sequence alignment files in various formats.

## Required Import

**Goal:** Load modules for reading, writing, and manipulating multiple sequence alignments.

**Approach:** Import AlignIO for file I/O and supporting classes for programmatic alignment construction.

```python
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
```

## Format Coverage Map

Three Python libraries cover the alignment-format space, with overlapping but non-identical support. Pick by what is actually required.

| Format | `Bio.AlignIO` | `Bio.Align` (modern) | `pyhmmer.easel` | Notes |
|--------|---------------|----------------------|-----------------|-------|
| Aligned FASTA | R/W | R/W | R/W | Most portable; loses annotations |
| Clustal | R/W | R/W | R | Clustal conservation marks NOT round-tripped |
| PHYLIP (interleaved/sequential/relaxed) | R/W | R/W | R | Strict 10-char names is silent footgun |
| Stockholm | R/W | R/W | R/W | Only format preserving GS/GR/GC/GF annotations |
| NEXUS | R/W | R/W | -- | MrBayes / PAUP* input |
| MAF (Multiple Alignment Format) | R/W | R/W | -- | UCSC whole-genome alignments |
| A2M / A3M | -- (use `'fasta'` parser then post-process) | -- | R/W | HMMER (a2m), HHsuite/ColabFold (a3m) |
| MSF (GCG) | R | -- | -- | GCG legacy |
| EMBOSS / Mauve XMFA / FASTA-m10 | R | partial | -- | One-way: read-only |

**Formats NOT in BioPython** (use dedicated tools):

| Format | Tool | Why |
|--------|------|-----|
| HAL | progressiveCactus, halTools | HDF5-backed multi-genome alignments at TB scale |
| chain / net | UCSC Kent tools (`liftOver`, `chainNet`) | Pairwise genome alignment |
| AXT | BLASTZ / lastz native | Pairwise alignment blocks |
| PSL | UCSC Kent tools (`pslPretty`, `blat`) | BLAT alignment summary |
| GFA / rGFA | `vg`, `odgi`, `pggb`, gfatools | Pangenome graph |
| GAF | `vg surject`, `vg call` | Graph alignment format (read-to-graph) |

Recommend `Bio.Align` (modern API) over `Bio.AlignIO` (legacy) for new code; it returns `Alignment` objects with built-in `.counts()` and `.substitutions` properties. For multi-gigabyte Stockholm databases such as Pfam-A.full, `pyhmmer.easel.MSAFile` streams record-by-record where `Bio.AlignIO.parse` works but at higher per-record cost.

## Reading Alignments

**"Read an alignment file"** → Parse an alignment file into an alignment object with sequences and metadata accessible.

**Goal:** Load alignment data from files in various formats (Clustal, PHYLIP, Stockholm, FASTA).

**Approach:** Use `AlignIO.read()` for single-alignment files or `AlignIO.parse()` for files containing multiple alignments.

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

**Goal:** Save alignment data to files in standard formats for downstream tools or archival.

**Approach:** Use `AlignIO.write()` with the target format specifier, supporting single or multiple alignments and file handles.

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

**"Convert alignment format"** → Transform an alignment file from one format to another (e.g., Clustal to PHYLIP).

**Goal:** Convert alignment files between formats for compatibility with different analysis tools.

**Approach:** Use `AlignIO.convert()` for direct one-step conversion, or read-modify-write for cases requiring intermediate manipulation.

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

**Goal:** Navigate and extract data from alignment objects including sequences, columns, and slices.

**Approach:** Use iteration, indexing, and column slicing on the alignment object.

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

**Goal:** Build an alignment object from sequences defined in code rather than read from a file.

**Approach:** Construct SeqRecord objects with gap characters and wrap them in a MultipleSeqAlignment.

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

## Format Selection for Downstream Tools

Choosing the output format depends on which downstream tool consumes the alignment:

| Downstream Tool | Required Format | BioPython Format String |
|----------------|-----------------|------------------------|
| RAxML-NG, IQ-TREE | PHYLIP (relaxed) | `'phylip-relaxed'` |
| MrBayes | NEXUS | `'nexus'` |
| PAUP* | NEXUS or PHYLIP | `'nexus'` or `'phylip'` |
| HMMER, Infernal | Stockholm | `'stockholm'` |
| Pfam/Rfam databases | Stockholm | `'stockholm'` |
| PAML/codeml | PHYLIP (sequential) | `'phylip-sequential'` |
| Most tools | FASTA | `'fasta'` |

### Annotation Preservation

Not all formats support annotations. Converting between formats can silently discard metadata:

| Format | Sequence Annotations | Column Annotations | Secondary Structure |
|--------|---------------------|-------------------|-------------------|
| Stockholm | Yes (GS/GR lines) | Yes (GC lines) | Yes (SS_cons) |
| NEXUS | Partial (SETS block) | Via CHARSET | No |
| Clustal | No (conservation marks not parsed) | No | No |
| PHYLIP | No | No | No |
| FASTA | No | No | No |

Converting Stockholm to FASTA or PHYLIP discards all annotations, secondary structure markup, and per-residue quality scores. If annotations matter, keep a Stockholm master copy.

## Format-Specific Notes

### PHYLIP Format Pitfalls

PHYLIP has two incompatible variants (interleaved vs sequential) and two name-length modes (strict vs relaxed). Confusing these causes silent data corruption.

**Strict PHYLIP** truncates sequence names to exactly 10 characters. This can silently merge distinct sequences whose names share a 10-character prefix (e.g., `Homo_sapiens_chr1` and `Homo_sapiens_chr2` both become `Homo_sapie`).

```python
# Strict PHYLIP (10-char names, interleaved) -- only for tools requiring it
alignment = AlignIO.read('file.phy', 'phylip')

# Sequential PHYLIP (10-char names, one sequence at a time) -- PAML/codeml
alignment = AlignIO.read('file.phy', 'phylip-sequential')

# Relaxed PHYLIP (no name limit) -- RAxML-NG, IQ-TREE (recommended default)
alignment = AlignIO.read('file.phy', 'phylip-relaxed')

# Always prefer phylip-relaxed for writing unless the downstream tool
# specifically requires strict format
AlignIO.write(alignment, 'output.phy', 'phylip-relaxed')
```

#### PHYLIP-Relaxed Dialect Mismatches Between Tree Tools

Biopython's `'phylip-relaxed'` writes a single space between name and sequence. RAxML-NG and IQ-TREE accept this; PhyML rejects sequence names containing colons or parentheses; PAML's codeml expects sequential format with name-truncation behaviour distinct from interleaved. Common silent failures:

| Symptom | Cause | Fix |
|---------|-------|-----|
| RAxML-NG: `terminating with uncaught exception ... bad alphabet` | Stop codons (`*`) in protein alignment | Replace `*` with `X` before writing |
| IQ-TREE: `not a valid PHYLIP file` | Sequence name contains `:` (NEXUS-tree-style refs) | Sanitize names: `re.sub(r'[():,]', '_', record.id)` |
| PhyML: silently truncated names | Names >100 chars | PhyML truncates without warning at 100 chars in current build |
| codeml: `cannot read sequences` | Used `phylip-relaxed` instead of `phylip-sequential` | codeml requires strict sequential |

Always verify by running the downstream tool's "validate input only" mode (e.g. `iqtree2 -s file.phy --check`) before committing to a long compute.

### MAF Block Coordinate Conventions

UCSC MAF (read via `AlignIO.parse(file, 'maf')`) returns blocks with per-row `annotations`:
- `start` (0-based; converts directly to BED but is off-by-one vs GFF)
- `size` (length on src strand)
- `strand` (`+` or `-`)
- `srcSize` (length of source chromosome)

For minus-strand rows, `start` is measured from the END of the source contig: the corresponding plus-strand start is `srcSize - start - size`. Without this conversion, lifting MAF to genome coordinates places minus-strand blocks at the wrong locus. Reference: UCSC MAF spec at genome.ucsc.edu/FAQ/FAQformat.html#format5.

```python
def maf_to_plus_strand_coords(row_anno):
    if row_anno['strand'] == '-':
        return row_anno['srcSize'] - row_anno['start'] - row_anno['size']
    return row_anno['start']
```

### Stockholm Format Annotations

Stockholm format (used by Pfam, Rfam, HMMER) supports four annotation line types:

| Line Prefix | Scope | Description | Example |
|-------------|-------|-------------|---------|
| `#=GF` | File | Alignment-level metadata (ID, accession, description) | `#=GF AC PF00001` |
| `#=GC` | Column | Per-column annotation (1 char per alignment column) | `#=GC SS_cons ..(((...)))..` |
| `#=GS` | Sequence | Per-sequence free text (organism, description) | `#=GS seq1 OS Homo sapiens` |
| `#=GR` | Residue | Per-residue annotation (1 char per residue) | `#=GR seq1 SS ..HHH..EEE..` |

Common GC annotations: `SS_cons` (consensus secondary structure), `RF` (reference coordinates), `seq_cons` (consensus sequence).

**WUSS notation** in RNA `#=GC SS_cons` lines uses nested bracket pairs (`<>`, `()`, `[]`, `{}`) for paired bases and characters like `_`, `-`, `,`, `:`, `.`, `~` for unpaired regions; pseudoknots use upper/lower-case letter pairs (`Aa`, `Bb`). Consult the Infernal user guide for the full character table before writing or parsing custom SS_cons strings.

```python
alignment = AlignIO.read('pfam.sto', 'stockholm')

for record in alignment:
    print(record.id, record.annotations)
    if 'secondary_structure' in record.letter_annotations:
        print(f'  SS: {record.letter_annotations["secondary_structure"]}')

ss_cons = alignment.column_annotations.get('secondary_structure')
```

**Round-trip caveat:** `AlignIO.write(alignment, 'out.fasta', 'fasta')` discards every Stockholm annotation silently. Re-reading and re-writing as Stockholm preserves GC/GR but dropped/added sequences invalidate the per-residue annotations -- regenerate annotations after edits.

**Pfam-style `name/start-end` identifier convention:** Pfam, Rfam, and Dfam Stockholm IDs (e.g. `Q9Y6Y0/45-198`) encode a 1-based inclusive region. Biopython does not split this; before passing to RAxML or IQ-TREE, parse the suffix into `record.annotations['start']` / `['end']` and strip from `record.id`, then restore it after.

### A2M / A3M Conventions

A2M (HMMER) and A3M (HHsuite, ColabFold) encode match vs insert columns by case (uppercase / `-` = match column, lowercase / `.` = insert column). A2M pads inserts across rows so it loads as a rectangular MSA; A3M does not, so convert with HHsuite `reformat.pl a3m a2m in.a3m out.a2m` (or `pyhmmer.easel.MSAFile(..., format='a2m')`) before parsing as a normal alignment.

**reformat.pl pitfall:** HHsuite's `reformat.pl a3m a2m` uses the FIRST sequence in the A3M as the match-state reference. ColabFold MSAs typically place the query first, which is the desired reference; merged or sorted A3Ms can have a non-query first sequence, producing match-state assignments that mis-align the query. Either (a) verify the first sequence is the query before reformatting, or (b) renormalise with `hhfilter -i in.a3m -o out.a3m -id 100 -qid 0 -cov 0` before running `reformat.pl`. A3M files emitted by `hhblits` always have the query first; A3M files concatenated from MSA databases do not.

```python
alignment = AlignIO.read('hhsearch.a2m', 'fasta')
match_only_seqs = [
    ''.join(c for c in str(r.seq) if c.isupper() or c == '-')
    for r in alignment
]
```

### Streaming Large Stockholm Databases

`Bio.AlignIO.read()` is in-memory; for Pfam-A.full (multi-gigabyte; ~22,000 family alignments in Pfam 37) or BFD (>2 TB), use `pyhmmer.easel.MSAFile` for streaming Stockholm or A3M.

```python
import pyhmmer

with pyhmmer.easel.MSAFile('Pfam-A.full', digital=True) as msa_file:
    for msa in msa_file:
        if msa.nseq < 50:
            continue
        weights = msa.compute_weights(method='pb')
        print(msa.name.decode(), msa.nseq, msa.alen, f'sum_w={sum(weights):.1f}')
```

`msa.compute_weights(method='pb')` computes Henikoff PB weights via the same Easel routine HMMER uses; the weights sum to the number of sequences (not Neff). For an Henikoff-style Neff estimate, see `msa-parsing/examples/neff.py`.

### Clustal Format
```python
# Clustal preserves conservation symbols in file but not when parsed
alignment = AlignIO.read('clustal.aln', 'clustal')
```

## Batch Processing Multiple Files

**Goal:** Convert a directory of alignment files from one format to another in bulk.

**Approach:** Glob for input files and iterate, reading each alignment and writing to the target format.

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

**Goal:** Use the modern Bio.Align module for alignment I/O with access to newer features like counts and substitutions.

**Approach:** Use `Align.read()`, `Align.parse()`, and `Align.write()` which return `Alignment` objects instead of `MultipleSeqAlignment`.

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

- alignment/multiple-alignment - Run MSA tools (MAFFT, MUSCLE5, ClustalOmega) to generate alignments
- alignment/pairwise-alignment - Create pairwise alignments with PairwiseAligner
- alignment/msa-parsing - Analyze alignment content and annotations
- alignment/msa-statistics - Calculate conservation and identity
- alignment/structural-alignment - Foldseek/TM-align outputs and Foldmason `result_aa.fa` / `result_3di.fa` MSAs (FASTA-loadable; the per-column LDDT report is HTML, not BioPython-parseable)
- alignment/alignment-trimming - Pre-format trimming with column-mapping retention
- sequence-io/format-conversion - Convert sequence (non-alignment) formats

## References

- Nawrocki EP, Eddy SR. 2013. Infernal 1.1: 100-fold faster RNA homology searches. Bioinf 29:2933-2935 (WUSS notation reference; see also the Infernal user guide).
- Larralde M et al. 2023. PyHMMER: a Python library binding to HMMER for efficient sequence analysis. Bioinf 39:btad214.
- Cock PJA et al. 2009. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinf 25:1422-1423.
- Mistry J et al. 2021. Pfam: the protein families database in 2021. NAR 49:D412-D419.
- Steinegger M et al. 2019. HH-suite3 for fast remote homology detection and deep protein annotation. BMC Bioinf 20:473.
- Mirdita M et al. 2022. ColabFold: making protein folding accessible to all. Nat Methods 19:679-682.
- Blanchette M et al. 2004. Aligning multiple genomic sequences with the threaded blockset aligner. Genome Res 14:708-715.

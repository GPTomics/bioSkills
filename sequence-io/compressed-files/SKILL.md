---
name: bio-compressed-files
description: Read, write, and index compressed sequence files (gzip, bzip2, xz, BGZF) with Biopython and bgzip/samtools. Use when working with .gz, .bz2, or .bgz sequence files, when random access into a compressed FASTA/FASTQ is needed, or when SeqIO.index/faidx/tabix rejects a plain .gz. Covers the BGZF-vs-gzip seekability asymmetry, the 'rt'-not-'rb' handle trap, virtual offsets, and gzip-to-BGZF conversion.
tool_type: mixed
primary_tool: Bio.bgzf
---

## Version Compatibility

Reference examples tested with: BioPython 1.83+, htslib/bgzip 1.19+, samtools 1.19+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Compressed Files

Read, write, and randomly access gzip, bzip2, xz, and BGZF compressed sequence files.

**"Read a compressed sequence file"** -> Open a decompression handle in TEXT mode, then parse with the standard SeqIO interface.
- gzip: `gzip.open(path, 'rt')` (Python stdlib)
- bzip2: `bz2.open(path, 'rt')` (Python stdlib)
- xz/LZMA: `lzma.open(path, 'rt')` (Python stdlib)
- BGZF: `bgzf.open(path, 'rt')` (BioPython) - BGZF input ONLY

**"Make a compressed file randomly accessible"** -> Re-compress as BGZF, then index. Only BGZF supports `SeqIO.index()`, `samtools faidx`, and `tabix` on compressed data.

## The Governing Principle: BGZF vs plain gzip is asymmetric

A BGZF (Blocked GNU Zip) file IS a valid gzip file - `gunzip` and `zcat` read it transparently. The reverse is FALSE: a plain `.gz` is NOT BGZF, so `faidx`/`tabix`/`SeqIO.index` reject it, and `bgzf.open` refuses to read it.

The reason is structural. BGZF is a series of concatenated gzip blocks, each <=64 KiB and independently decodable, so any record can be reached by seeking to its block. Plain gzip is one continuous DEFLATE stream with no block boundaries: reaching byte N means decompressing every byte before it (O(n)). Random access therefore REQUIRES BGZF; on a plain `.gz`, `SeqIO.index()` would re-decompress huge prefixes on every lookup, which is why Biopython forbids it outright (it raises rather than running slowly).

Consequences the agent must respect:
- Reading sequentially: gzip, bzip2, xz, and BGZF all work via the matching `*.open(path, 'rt')` handle.
- Random access / indexing: BGZF only. Convert plain gzip to BGZF first.
- `bgzf.open()` reads BGZF input only. Pointing it at a plain `.gz` raises `ValueError: A BGZF block should start with b'\x1f\x8b\x08\x04'...`. To read plain gzip use `gzip.open()`.
- BAM and tabix-indexed files use BGZF natively; bzip2/xz are archive-only (no seekable index).

## Required Imports

```python
import gzip
import bz2
import lzma
from Bio import SeqIO
from Bio import bgzf
```

## Reading Compressed Files

**Goal:** Parse sequence records from a compressed file without decompressing to disk.

**Approach:** Open a decompression handle in TEXT mode (`'rt'`), then pass the handle to `SeqIO.parse()`. The parser is format-agnostic about the underlying compression.

```python
with gzip.open('reads.fastq.gz', 'rt') as handle:
    for record in SeqIO.parse(handle, 'fastq'):
        print(record.id, len(record.seq))
```

Swap `gzip.open` for `bz2.open` (`.bz2`), `lzma.open` (`.xz`), or `bgzf.open` (`.bgz`) - the parse loop is identical.

### The 'rt' vs 'rb' trap

`SeqIO.parse()` in Python 3 needs a TEXT handle that yields `str`. A binary `'rb'` handle yields `bytes` and raises `TypeError: a bytes-like object is required` (or a decode error). Always use `'rt'` for reading and `'wt'` for writing through `SeqIO`. The low-level `SimpleFastaParser`/`FastqGeneralIterator` also require text handles.

## Writing Compressed Files

**Goal:** Save records straight to a compressed file with no intermediate plain copy.

**Approach:** Open a compression handle in TEXT mode (`'wt'`), then pass it to `SeqIO.write()`.

```python
with gzip.open('output.fasta.gz', 'wt') as handle:
    SeqIO.write(records, handle, 'fasta')
```

For an indexable result write BGZF instead:

```python
with bgzf.open('output.fasta.bgz', 'wt') as handle:
    SeqIO.write(records, handle, 'fasta')
```

`BgzfWriter.close()` (and the `with` block exit) automatically appends the 28-byte empty-block EOF marker that htslib tools check for; let the context manager close the handle.

## Random Access: index a BGZF file

**Goal:** Pull individual records by id from a large compressed file without a linear scan.

**Approach:** Compress as BGZF, then build a Biopython offset index. `SeqIO.index()` keeps virtual offsets in RAM; `SeqIO.index_db()` stores them in an on-disk SQLite index that persists across sessions and spans multiple files.

```python
records = SeqIO.index('sequences.fasta.bgz', 'fasta')
target = records['gene_042'].seq
records.close()

# Persistent, multi-file, scales beyond RAM:
db = SeqIO.index_db('idx.sqlite', ['a.fasta.bgz', 'b.fasta.bgz'], 'fasta')
```

`SeqIO.index()` on a plain `.gz` raises `ValueError: Gzipped files are not suitable for indexing, please use BGZF (blocked gzip format) instead.` Convert first (below).

## Convert plain gzip to BGZF

**"Convert gzip to an indexable format"** -> Decompress the gzip stream and re-compress it as BGZF.

CLI (fastest, htslib-native):

```bash
# Either decompress then bgzip in place...
gzip -d sequences.fasta.gz && bgzip sequences.fasta      # -> sequences.fasta.gz (now BGZF)
# ...or stream without touching disk:
zcat sequences.fasta.gz | bgzip -@ 4 > sequences.fasta.bgz

# Index a BGZF FASTA for region extraction:
samtools faidx sequences.fasta.bgz                        # writes BOTH .fai and .gzi
samtools faidx sequences.fasta.bgz gene_042:1-200
```

`samtools faidx` on a BGZF FASTA writes TWO index files: `.fai` (record offsets in uncompressed coordinates) AND `.gzi` (the compressed-to-uncompressed block map). Deleting `.gzi` breaks region extraction even though `.fai` survives. Note that bgzip keeps the `.gz` extension, so a `.gz` may be EITHER plain gzip or BGZF - check with `bgzip -t file.gz` (tests for a valid BGZF stream) rather than trusting the suffix.

Pure-Python equivalent (no external tools):

```python
import gzip
from Bio import SeqIO, bgzf

with gzip.open('input.fasta.gz', 'rt') as in_handle:
    with bgzf.open('output.fasta.bgz', 'wt') as out_handle:
        SeqIO.write(SeqIO.parse(in_handle, 'fasta'), out_handle, 'fasta')
```

## Bio.bgzf API and virtual offsets

`Bio.bgzf` exports `open`, `BgzfReader`, `BgzfWriter`, `make_virtual_offset`, `split_virtual_offset`.

A virtual offset packs two coordinates into one 64-bit integer: `voffset = coffset << 16 | uoffset`, where `coffset` is the byte position of the block start in the compressed file (top 48 bits) and `uoffset` is the offset within that block's decompressed data (low 16 bits - 16 bits suffices because a block holds at most 64 KiB).

```python
vo = bgzf.make_virtual_offset(100, 7)        # 6553607
coffset, uoffset = bgzf.split_virtual_offset(vo)   # (100, 7)

with bgzf.open('sequences.fasta.bgz', 'rt') as handle:
    handle.readline()
    saved = handle.tell()        # a VIRTUAL offset, not a byte position
    handle.seek(saved)           # jumps back to the same record
```

Critical caveat: virtual offsets may be COMPARED for ordering but NEVER SUBTRACTED to get a byte length - they live in two coordinate spaces (compressed position and within-block position), so `vo2 - vo1` is meaningless. `BgzfReader.tell()` returns a virtual offset; `seek()` consumes one. Text mode forces `latin1` and does no newline translation.

## Compression Format Comparison

| Format | Extension | Random access | Speed | Ratio | Stdlib handle |
|--------|-----------|---------------|-------|-------|---------------|
| gzip | `.gz` | No (O(n) seek) | Fast | Good | `gzip.open` |
| BGZF | `.bgz` / `.gz` | **Yes (block-seekable)** | Fast, threadable | Good | `bgzf.open` (BioPython) |
| bzip2 | `.bz2` | No | Slow | Better | `bz2.open` |
| xz / LZMA | `.xz` | No | Slowest | Best | `lzma.open` |

## When to Use Each Format

| Use case | Format | Why |
|----------|--------|-----|
| Sequential read/write, sharing | gzip | Universal, fast, every tool reads it |
| Need `faidx`/`tabix`/`SeqIO.index` | **BGZF** | Only seekable compressed format |
| BAM, tabix-indexed VCF/GFF/BED | BGZF | Required natively |
| Cold archive, max shrink, no random access | xz then bzip2 | Highest ratios, slowest |
| Random access into an existing plain `.gz` without re-bgzipping | pyfastx | Adds a seek-point index over the gzip stream |

`pyfastx` is the exception that gives random access into a PLAIN gzip FASTA/FASTQ: it builds a seek-point index (via `zran` from indexed_gzip) plus a SQLite `.fxi`/`.fqi` index alongside the file, a different strategy from faidx (which requires the stream itself to be BGZF). Use it when re-compressing a large gzipped genome to BGZF is not an option.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| `TypeError: a bytes-like object is required` | Handle opened `'rb'` instead of `'rt'` | Open compressed handles with `'rt'`/`'wt'` for SeqIO |
| `ValueError: Gzipped files are not suitable for indexing, please use BGZF...` | `SeqIO.index()` on a plain `.gz` | Re-compress as BGZF (`zcat ... | bgzip`), then index |
| `ValueError: A BGZF block should start with b'\x1f\x8b\x08\x04'...` | `bgzf.open()` pointed at a plain gzip file | Read plain gzip with `gzip.open()`; reserve `bgzf.open` for BGZF |
| `[bgzf] file ... not BGZF` / `not compressed with bgzip` (htslib) | faidx/tabix given a plain `.gz` | Convert to BGZF first |
| `[faidx] Failed to read ... / could not load .gzi` | `.gzi` deleted next to a BGZF FASTA | Re-run `samtools faidx` to regenerate `.fai` + `.gzi` |
| `gzip.BadGzipFile` / `OSError: Not a gzipped file` | File is not gzip (wrong suffix / corrupt) | Verify with `bgzip -t` or `file`; match handle to real format |
| `UnicodeDecodeError` | Non-UTF8 bytes in a text handle | `gzip.open(path, 'rt', encoding='latin-1')` |

## Related Skills

- read-sequences - parse vs index vs index_db trade-offs for compressed handles
- write-sequences - write records through a compression handle
- batch-processing - stream many compressed files without loading them into RAM
- filter-sequences - keep R1/R2 in sync when filtering gzipped paired reads
- alignment-files/sam-bam-basics - BAM is BGZF natively; samtools manages the compression

## References

- Li H, Handsaker B, Wysoker A, et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009;25(16):2078-2079. (Defines BGZF in the SAM/BAM specification.)
- Bonfield JK. CRAM 3.1: advances in the CRAM file format. Bioinformatics. 2022;38(6):1497. (Per-column block compression beyond BGZF.)
- Du L, Liu Q, Fan Z, et al. Pyfastx: a robust Python package for fast random access to sequences from plain and gzipped FASTA/Q files. Briefings in Bioinformatics. 2021;22(4):bbaa368.

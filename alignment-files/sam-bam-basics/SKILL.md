---
name: bio-sam-bam-basics
description: View, convert, and understand SAM/BAM/CRAM alignment files using samtools and pysam. Use when inspecting alignments, converting between formats, or understanding alignment file structure.
tool_type: cli
primary_tool: samtools
---

## Version Compatibility

Reference examples tested with: pysam 0.22+, samtools 1.19+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# SAM/BAM/CRAM Basics

**"Read a BAM file"** → Open a binary alignment file and iterate over aligned reads with their mapping coordinates, flags, and quality scores.
- Python: `pysam.AlignmentFile()` (pysam)
- CLI: `samtools view` (samtools)
- R: `scanBam()` (Rsamtools)

View and convert alignment files using samtools and pysam.

## Format Overview

| Format | Description | Use Case |
|--------|-------------|----------|
| SAM | Text format, human-readable | Debugging, small files |
| BAM | Binary compressed SAM | Standard storage format |
| CRAM | Reference-based compression | Long-term archival, smaller than BAM |

## SAM Format Structure

```
@HD VN:1.6 SO:coordinate
@SQ SN:chr1 LN:248956422
@RG ID:sample1 SM:sample1
@PG ID:bwa PN:bwa VN:0.7.17
read1  0   chr1  100  60  50M  *  0  0  ACGT...  FFFF...  NM:i:0
```

Header lines start with `@`:
- `@HD` - Header metadata (version, sort order)
- `@SQ` - Reference sequence dictionary
- `@RG` - Read group information
- `@PG` - Program used to create file

Alignment fields (tab-separated):
1. QNAME - Read name
2. FLAG - Bitwise flag
3. RNAME - Reference name
4. POS - 1-based position
5. MAPQ - Mapping quality
6. CIGAR - Alignment description
7. RNEXT - Mate reference
8. PNEXT - Mate position
9. TLEN - Template length
10. SEQ - Read sequence
11. QUAL - Base qualities
12. Optional tags (NM:i:0, MD:Z:50, etc.)

## samtools view

### View BAM as SAM
```bash
samtools view input.bam | head
```

### View with Header
```bash
samtools view -h input.bam | head -100
```

### View Header Only
```bash
samtools view -H input.bam
```

### View Specific Region
```bash
samtools view input.bam chr1:1000-2000
```

### Count Alignments
```bash
samtools view -c input.bam
```

## Format Conversion

**Goal:** Convert between SAM (text), BAM (binary), and CRAM (reference-compressed) alignment formats.

**Approach:** Use `samtools view` with format flags (`-b` for BAM, `-C` for CRAM, `-h` for SAM with header). CRAM requires a reference FASTA with `-T`.

### BAM to SAM
```bash
samtools view -h -o output.sam input.bam
```

### SAM to BAM
```bash
samtools view -b -o output.bam input.sam
```

### BAM to CRAM
```bash
samtools view -C -T reference.fa -o output.cram input.bam
```

### CRAM to BAM
```bash
samtools view -b -T reference.fa -o output.bam input.cram
```

### Pipe Conversion
```bash
samtools view -b input.sam > output.bam
```

## Common Flags

| Flag | Decimal | Meaning |
|------|---------|---------|
| 0x1 | 1 | Paired |
| 0x2 | 2 | Proper pair |
| 0x4 | 4 | Unmapped |
| 0x8 | 8 | Mate unmapped |
| 0x10 | 16 | Reverse strand |
| 0x20 | 32 | Mate reverse strand |
| 0x40 | 64 | First in pair |
| 0x80 | 128 | Second in pair |
| 0x100 | 256 | Secondary alignment |
| 0x200 | 512 | Failed QC |
| 0x400 | 1024 | PCR duplicate |
| 0x800 | 2048 | Supplementary |

### Decode Flags (Bidirectional)
```bash
# Number to mnemonics
samtools flags 147
# 0x93 147 PAIRED,PROPER_PAIR,REVERSE,READ2

# Mnemonics to number
samtools flags PAIRED,PROPER_PAIR,REVERSE,READ2   # 147
```

### Secondary vs Supplementary (Different Semantics)

Two different concepts that are routinely conflated:

| Bit | Name | Meaning | Filter implication |
|-----|------|---------|--------------------|
| 0x100 (256) | Secondary | An alternative candidate alignment for the same read; not the primary location | `-F 256` is correct for SNV/indel calling on short reads |
| 0x800 (2048) | Supplementary | A piece of a chimeric/split alignment (the read is split across loci) | Carries SA:Z tag; **required** by SV callers (Manta, Sniffles, cuteSV, GRIDSS, Delly) |

`-F 2304` removes both. Strip supplementary only when downstream is small-variant calling; keep supplementary for SV calling, fusion detection, or any analysis that follows split-reads.

## MAPQ Is Not Portable Across Aligners

`samtools view -q 30` does different things depending on what produced the BAM. MAPQ is an aligner-specific scale, not a universal probability:

| Aligner | MAPQ scale | "Unique" sentinel | Common gotcha |
|---------|-----------|-------------------|----------------|
| BWA-MEM / BWA-MEM2 | 0-60 | 60 | `-q 30` is sensible "high confidence" |
| minimap2 (DNA / pbmm2) | 0-60 | 60 | Spec-compliant |
| HISAT2 | 0-60 | 60 | Spec-compliant |
| Bowtie2 | 0-42 | 42 (rare) | `-q 60` drops everything; `-q 23` is the established 99% threshold |
| STAR | 0, 1, 2, 3, 255 | **255 = uniquely mapped (sentinel, not a quality)** | `-q 255` for "unique only"; `-q 30` accidentally keeps unique only too |
| DRAGEN | 0 to `--mapq-max` (often ~250) | varies | `-q 30` still meaningful; distribution shape differs |
| Cell Ranger / STARsolo | inherits STAR | 255 | Same trap as STAR |

Verify the actual scale of any unfamiliar BAM:
```bash
samtools view input.bam | awk '{print $5}' | sort -un | head
samtools view -H input.bam | grep '^@PG' | head -1   # which aligner produced this BAM
```

## 0-Based vs 1-Based Coordinates (Footgun)

| Context | Coordinate system |
|---------|-------------------|
| SAM text POS | 1-based, inclusive |
| `samtools view chr1:100-200` | 1-based, closed interval |
| `samtools faidx chr1:100-200` | 1-based, closed interval |
| BAM binary internal | 0-based, half-open |
| `pysam read.reference_start` | 0-based |
| `bam.fetch('chr1', 100, 200)` | 0-based, half-open |
| BED files | 0-based, half-open |
| VCF | 1-based |
| GFF/GTF | 1-based, inclusive |

`samtools view bam chr1:100-200` and `bam.fetch('chr1', 100, 200)` return different read sets at boundaries.

## CIGAR Operations

| Op | Description |
|----|-------------|
| M | Alignment match (can be mismatch) |
| I | Insertion to reference |
| D | Deletion from reference |
| N | Skipped region (introns in RNA-seq; do NOT count as covered bases) |
| S | Soft clipping (sequence in SEQ but not aligned) |
| H | Hard clipping (sequence not in SEQ) |
| = | Sequence match (explicit) |
| X | Sequence mismatch (explicit) |
| P | Padding (rare; multiple-sequence-alignment context) |

Example: `50M2I30M` = 50 bases match, 2 base insertion, 30 bases match

CIGAR `M` is overloaded -- it is the union of `=` and `X`. Some aligners (minimap2, BWA with `-Y`) emit `=`/`X` directly; bcftools / Picard often need `M` and rebuild MD/NM with `samtools calmd`. `N` operations break naive coverage calculations: a 1000 bp RNA-seq read with one 50 kb intron does not cover 50 kb. Distinguish soft-clip (`S`, bases retained) from hard-clip (`H`, bases discarded -- irreversible).

## Context-Specific Tags

Beyond the standard fields, downstream tools depend on optional tags whose presence depends on aligner and assay. Inspect with `samtools view input.bam | head -1 | tr '\t' '\n'` or pysam `read.get_tag('XX')`.

| Tag | Set by | Meaning | Required by |
|-----|--------|---------|-------------|
| NM:i | bwa, samtools calmd | Edit distance to reference | mapDamage, many filters |
| MD:Z | bwa, samtools calmd | Mismatch positions (text) | bcftools mpileup BAQ, IGV mismatch coloring |
| MC:Z | samtools fixmate -m | Mate CIGAR | samtools markdup |
| MS:i | samtools fixmate -m | Mate score | samtools markdup |
| RG:Z | aligner from -R | Read group ID | GATK BQSR, MarkDuplicates LB lookup |
| SA:Z | All split-read aligners | Comma-list of supplementary coords | Sniffles, Manta, cuteSV, GRIDSS, Delly |
| NH:i | STAR, HISAT2 | Number of reported hits | featureCounts multimapper handling, Salmon |
| HI:i | STAR | Hit index (0-based among NH) | RSEM |
| XS:A | STAR, HISAT2, minimap2 -ax splice | Strand inferred from splice motif | StringTie, Cufflinks |
| CB:Z | Cell Ranger, STARsolo | Corrected cell barcode | scRNA quantification |
| UB:Z | Cell Ranger, STARsolo | Corrected UMI | UMI-aware dedup |
| RX:Z | fgbio AnnotateBamWithUmis | Raw UMI (bulk) | fgbio GroupReadsByUmi |
| MI:Z | fgbio CallMolecularConsensusReads | Molecular identifier (consensus) | Duplex calling |
| cs:Z | minimap2 --cs | Compact CIGAR-with-bases | paftools, SV tools |

Missing tags fail in two modes: silently wrong (featureCounts ignoring multimappers without NH; markdup marking nothing without MC/MS) or loudly (consensus tools rejecting input without MD).

## Provenance: @PG Chain

The `@PG` lines record every tool that touched the BAM, linked through `PP` (previous program) tags. This is the audit trail.

```bash
samtools view -H input.bam | grep '^@PG'
```

A clean germline pipeline:
```
@PG ID:bwa-mem PN:bwa VN:0.7.17
@PG ID:samtools.1 PN:samtools VN:1.20 PP:bwa-mem CL:samtools sort
@PG ID:samtools.2 PN:samtools VN:1.20 PP:samtools.1 CL:samtools fixmate
@PG ID:samtools.3 PN:samtools VN:1.20 PP:samtools.2 CL:samtools markdup
```

A broken/missing chain (no PP, unknown tools, gaps) means the BAM cannot be reliably reproduced. Production pipelines often reject inputs without a complete chain.

## CRAM Reference Resolution (Critical)

CRAM stores reads relative to a reference; without it, the file is unreadable. htslib resolves the reference in this order:

1. Command-line `-T ref.fa` / `--reference`
2. `REF_CACHE` env var (local MD5-named cache)
3. `REF_PATH` env var (colon-separated; can include URLs)
4. `UR:` URL in SAM `@SQ` header
5. Last resort: EBI ENA download via MD5 in `M5:` tag (fails on offline HPC)

On HPC nodes without internet, populate a local cache once:
```bash
mkdir -p $HOME/cram_cache
seq_cache_populate.pl -root $HOME/cram_cache reference.fa
export REF_CACHE=$HOME/cram_cache/%2s/%2s/%s
export REF_PATH=$REF_CACHE   # disables ENA fallback

samtools quickcheck -v file.cram   # header + EOF only
samtools view -c file.cram          # forces full decode; proves reference reachable
```

CRAM operations can be irreversibly lossy: `--output-fmt-option=archive=1` enables 8-bin Illumina quality binning (~30-50% additional size reduction; benign for >=30x germline WGS, harmful for low-coverage / somatic / forensic / archival). Convert against the *exact* reference the BAM was aligned to (matched by `@SQ M5:`); a different reference silently corrupts bases on read-back.

## pysam Python Alternative

**Goal:** Read and manipulate alignment data programmatically in Python.

**Approach:** Use `pysam.AlignmentFile` to open BAM/CRAM files, iterate over reads, and access properties like coordinates, flags, CIGAR, and tags.

### Open and Iterate
```python
import pysam

with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for read in bam:
        print(f'{read.query_name}\t{read.reference_name}:{read.reference_start}')
```

### Access Header
```python
with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for sq in bam.header['SQ']:
        print(f'{sq["SN"]}: {sq["LN"]} bp')
```

### Read Alignment Properties
```python
with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for read in bam:
        print(f'Name: {read.query_name}')
        print(f'Flag: {read.flag}')
        print(f'Chrom: {read.reference_name}')
        print(f'Pos: {read.reference_start}')  # 0-based
        print(f'MAPQ: {read.mapping_quality}')
        print(f'CIGAR: {read.cigarstring}')
        print(f'Seq: {read.query_sequence}')
        print(f'Qual: {read.query_qualities}')
        break
```

### Check Flag Properties
```python
with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for read in bam:
        if read.is_paired and read.is_proper_pair:
            if read.is_reverse:
                strand = '-'
            else:
                strand = '+'
            print(f'{read.query_name} on {strand} strand')
```

### Fetch Region
```python
with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for read in bam.fetch('chr1', 1000, 2000):
        print(read.query_name)
```

### Convert BAM to SAM
```python
with pysam.AlignmentFile('input.bam', 'rb') as infile:
    with pysam.AlignmentFile('output.sam', 'w', header=infile.header) as outfile:
        for read in infile:
            outfile.write(read)
```

### Convert to CRAM
```python
with pysam.AlignmentFile('input.bam', 'rb') as infile:
    with pysam.AlignmentFile('output.cram', 'wc', reference_filename='reference.fa', header=infile.header) as outfile:
        for read in infile:
            outfile.write(read)
```

## Quick Reference

| Task | samtools | pysam |
|------|----------|-------|
| View BAM | `samtools view file.bam` | `AlignmentFile('file.bam', 'rb')` |
| View header | `samtools view -H file.bam` | `bam.header` |
| Count reads | `samtools view -c file.bam` | `sum(1 for _ in bam)` |
| Get region | `samtools view file.bam chr1:1-1000` | `bam.fetch('chr1', 0, 1000)` |
| BAM to SAM | `samtools view -h -o out.sam in.bam` | Open with 'w' mode |
| SAM to BAM | `samtools view -b -o out.bam in.sam` | Open with 'wb' mode |
| BAM to CRAM | `samtools view -C -T ref.fa -o out.cram in.bam` | Open with 'wc' mode |

## Related Skills

- alignment-indexing - Create indices for random access (required for fetch/region queries)
- alignment-sorting - Sort alignments by coordinate or name
- alignment-filtering - Filter alignments by flags, quality, regions
- alignment-validation - Sequence dictionary cross-validation (M5 checksums)
- bam-statistics - Generate statistics from alignment files
- reference-operations - REF_PATH/REF_CACHE setup for CRAM
- sequence-io/read-sequences - Parse FASTA/FASTQ input files

# SAM/BAM/CRAM Basics Usage Guide

This guide covers viewing, converting, and understanding alignment files.

## Prerequisites

- samtools installed (`conda install -c bioconda samtools` or `brew install samtools`)
- pysam installed (`pip install pysam`)

## Understanding SAM Format

SAM (Sequence Alignment/Map) is the standard text format for storing aligned sequencing reads. BAM is the compressed binary equivalent, and CRAM provides even better compression using reference-based encoding.

### File Structure

A SAM file has two sections:

1. **Header** - Lines starting with `@` containing metadata
2. **Alignments** - Tab-separated alignment records

### Header Types

```
@HD VN:1.6 SO:coordinate           # File metadata
@SQ SN:chr1 LN:248956422           # Reference sequences
@RG ID:sample1 SM:sample1 PL:ILLUMINA  # Read groups
@PG ID:bwa PN:bwa VN:0.7.17        # Programs used
```

### Alignment Fields

Each alignment line has 11 mandatory fields plus optional tags:

| Column | Name | Description |
|--------|------|-------------|
| 1 | QNAME | Query (read) name |
| 2 | FLAG | Bitwise flags encoding read properties |
| 3 | RNAME | Reference sequence name |
| 4 | POS | 1-based leftmost mapping position |
| 5 | MAPQ | Mapping quality (Phred-scaled) |
| 6 | CIGAR | Alignment string |
| 7 | RNEXT | Reference name of mate/next read |
| 8 | PNEXT | Position of mate/next read |
| 9 | TLEN | Observed template length |
| 10 | SEQ | Read sequence |
| 11 | QUAL | Base quality string |

## Working with samtools

### Viewing Files

```bash
# View first 10 alignments
samtools view input.bam | head

# View with header
samtools view -h input.bam | head -50

# View header only
samtools view -H input.bam

# View specific chromosome
samtools view input.bam chr1

# View specific region (requires index)
samtools view input.bam chr1:1000000-2000000

# Count total alignments
samtools view -c input.bam

# Count alignments in region
samtools view -c input.bam chr1:1000000-2000000
```

### Format Conversion

```bash
# SAM to BAM (compression)
samtools view -b -o output.bam input.sam

# BAM to SAM (for inspection/debugging)
samtools view -h -o output.sam input.bam

# BAM to CRAM (maximum compression, needs reference)
samtools view -C -T reference.fa -o output.cram input.bam

# CRAM to BAM
samtools view -b -T reference.fa -o output.bam input.cram
```

### Understanding FLAGS

The FLAG field is a bitwise combination of properties:

```bash
# Decode a FLAG value
samtools flags 99
# Output: 0x63 99 PAIRED,PROPER_PAIR,MREVERSE,READ1

samtools flags 147
# Output: 0x93 147 PAIRED,PROPER_PAIR,REVERSE,READ2
```

Common FLAG combinations:
- 99: First read, properly paired, mate on reverse strand
- 147: Second read, properly paired, on reverse strand
- 4: Unmapped read
- 256: Secondary alignment

## Working with pysam

### Basic Reading

```python
import pysam

# Open BAM file
bam = pysam.AlignmentFile('input.bam', 'rb')

# Iterate over all reads
for read in bam:
    print(read.query_name, read.reference_name, read.reference_start)

bam.close()
```

### Using Context Manager

```python
import pysam

with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for read in bam:
        # Process read
        pass
```

### Accessing Read Properties

```python
import pysam

with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for read in bam:
        # Basic properties
        name = read.query_name
        flag = read.flag
        chrom = read.reference_name
        pos = read.reference_start  # 0-based!
        mapq = read.mapping_quality
        cigar = read.cigarstring
        seq = read.query_sequence
        qual = read.query_qualities  # As array of integers

        # Flag checks
        is_paired = read.is_paired
        is_proper = read.is_proper_pair
        is_unmapped = read.is_unmapped
        is_reverse = read.is_reverse
        is_read1 = read.is_read1
        is_read2 = read.is_read2
        is_secondary = read.is_secondary
        is_supplementary = read.is_supplementary
        is_duplicate = read.is_duplicate

        # Optional tags
        nm = read.get_tag('NM') if read.has_tag('NM') else None

        break  # Just show first read
```

### Fetching Regions

```python
import pysam

with pysam.AlignmentFile('input.bam', 'rb') as bam:
    # Fetch requires index (input.bam.bai)
    for read in bam.fetch('chr1', 1000000, 2000000):
        print(read.query_name)
```

### Writing Files

```python
import pysam

# Copy BAM to new file
with pysam.AlignmentFile('input.bam', 'rb') as infile:
    with pysam.AlignmentFile('output.bam', 'wb', header=infile.header) as outfile:
        for read in infile:
            outfile.write(read)
```

### File Mode Strings

| Mode | Description |
|------|-------------|
| `rb` | Read BAM |
| `r` | Read SAM |
| `rc` | Read CRAM |
| `wb` | Write BAM |
| `w` | Write SAM |
| `wc` | Write CRAM |
| `wbu` | Write uncompressed BAM |

## CIGAR String Interpretation

CIGAR describes how the read aligns to the reference:

| Operation | Consumes Query | Consumes Reference | Description |
|-----------|----------------|-------------------|-------------|
| M | Yes | Yes | Alignment match/mismatch |
| I | Yes | No | Insertion |
| D | No | Yes | Deletion |
| N | No | Yes | Skipped region (RNA intron) |
| S | Yes | No | Soft clipping |
| H | No | No | Hard clipping |
| = | Yes | Yes | Sequence match |
| X | Yes | Yes | Sequence mismatch |

Example: `10M2I30M5D20M`
- 10 bases match
- 2 bases inserted (in read, not in reference)
- 30 bases match
- 5 bases deleted (in reference, not in read)
- 20 bases match

## Common Workflows

### Inspect a BAM File

```bash
# Check if sorted and indexed
samtools view -H input.bam | grep "^@HD"

# Count reads
samtools view -c input.bam

# View first few reads
samtools view input.bam | head

# Check reference sequences
samtools view -H input.bam | grep "^@SQ" | head
```

### Convert Aligner Output

```bash
# BWA outputs SAM, convert to sorted BAM
bwa mem reference.fa reads.fq | samtools sort -o aligned.bam

# Index the result
samtools index aligned.bam
```

### Prepare for IGV

```bash
# IGV needs sorted, indexed BAM
samtools sort -o sorted.bam input.bam
samtools index sorted.bam
```

## Troubleshooting

| Issue | Cause | Solution |
|-------|-------|----------|
| `[E::idx_find_and_load] Could not retrieve index file` | Missing index | Run `samtools index file.bam` |
| `[E::hts_open_format] Failed to open file` | File doesn't exist or wrong permissions | Check path and permissions |
| Truncated output | Pipe closed early | Use `samtools view -h` to include header |
| CRAM errors | Missing reference | Provide `-T reference.fa` |

## See Also

- [SAM Format Specification](https://samtools.github.io/hts-specs/SAMv1.pdf)
- [samtools documentation](http://www.htslib.org/doc/samtools.html)
- [pysam documentation](https://pysam.readthedocs.io/)

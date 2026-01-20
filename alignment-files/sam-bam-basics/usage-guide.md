# SAM/BAM/CRAM Basics - Usage Guide

## Overview
View, convert, and understand alignment files in SAM, BAM, and CRAM formats using samtools and pysam.

## Prerequisites
```bash
# samtools
conda install -c bioconda samtools

# pysam
pip install pysam
```

## Quick Start
Tell your AI agent what you want to do:
- "View the first 10 reads in my BAM file"
- "Convert my SAM file to BAM format"
- "Count the total number of reads in sample.bam"
- "Extract reads from chromosome 1"

## Example Prompts

### Viewing Files
> "Show me the header of my BAM file"

> "View the first 20 alignments with the header included"

> "Count how many reads are in my BAM file"

### Format Conversion
> "Convert my SAM file to compressed BAM"

> "Convert my BAM to CRAM format using the reference genome"

> "Convert CRAM back to BAM for compatibility with older tools"

### Region Extraction
> "Extract all reads from chr1:1000000-2000000"

> "Get reads from multiple regions: chr1:1000-2000 and chr2:3000-4000"

> "Count reads in a specific genomic region"

### Understanding Alignments
> "Explain what FLAG 99 means in my SAM file"

> "Parse the CIGAR string 10M2I30M5D20M"

> "Show me the mapping quality distribution of my reads"

## What the Agent Will Do

1. Open the alignment file with appropriate mode (SAM/BAM/CRAM)
2. Parse the header to understand reference sequences and metadata
3. Iterate through alignments or fetch specific regions
4. Extract read properties (name, flag, position, quality, sequence)
5. Apply any requested filters or transformations
6. Output results in the requested format

## Understanding SAM Format

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

### Common FLAGS
| FLAG | Meaning |
|------|---------|
| 99 | First read, properly paired, mate on reverse strand |
| 147 | Second read, properly paired, on reverse strand |
| 4 | Unmapped read |
| 256 | Secondary alignment |
| 2048 | Supplementary alignment |

### CIGAR Operations
| Op | Description |
|----|-------------|
| M | Alignment match/mismatch |
| I | Insertion to reference |
| D | Deletion from reference |
| N | Skipped region (RNA intron) |
| S | Soft clipping |
| H | Hard clipping |

## Common Commands

### Viewing
```bash
samtools view input.bam | head              # First 10 alignments
samtools view -h input.bam | head -50       # With header
samtools view -H input.bam                  # Header only
samtools view input.bam chr1                # Specific chromosome
samtools view input.bam chr1:1000000-2000000  # Specific region
samtools view -c input.bam                  # Count alignments
```

### Conversion
```bash
samtools view -b -o output.bam input.sam              # SAM to BAM
samtools view -h -o output.sam input.bam              # BAM to SAM
samtools view -C -T reference.fa -o output.cram input.bam  # BAM to CRAM
samtools view -b -T reference.fa -o output.bam input.cram  # CRAM to BAM
```

### Decode FLAGS
```bash
samtools flags 99   # 0x63 99 PAIRED,PROPER_PAIR,MREVERSE,READ1
samtools flags 147  # 0x93 147 PAIRED,PROPER_PAIR,REVERSE,READ2
```

## Python with pysam

### Reading BAM
```python
import pysam

with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for read in bam:
        print(read.query_name, read.reference_name, read.reference_start)
```

### Accessing Read Properties
```python
import pysam

with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for read in bam:
        name = read.query_name
        flag = read.flag
        chrom = read.reference_name
        pos = read.reference_start  # 0-based
        mapq = read.mapping_quality
        cigar = read.cigarstring
        seq = read.query_sequence
        qual = read.query_qualities
        break
```

### Fetching Regions
```python
import pysam

with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for read in bam.fetch('chr1', 1000000, 2000000):
        print(read.query_name)
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

## Troubleshooting

### Missing Index Error
```bash
# Error: Could not retrieve index file
samtools index file.bam
```

### CRAM Requires Reference
```bash
# Always provide reference for CRAM
samtools view -T reference.fa input.cram
```

### Chromosome Name Mismatch
```bash
# Check chromosome names in file
samtools view -H input.bam | grep "^@SQ"
# Use matching names: "chr1" not "1" if file uses "chr1"
```

## Tips
- Always check if BAM is sorted and indexed before region queries
- Use BAM for most operations, CRAM only when storage is critical
- pysam uses 0-based coordinates, samtools uses 1-based
- Include `-h` flag when piping to preserve header information
- Use context managers with pysam to ensure proper file closing
- The `-@` flag enables multi-threading for faster processing

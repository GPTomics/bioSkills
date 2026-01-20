# Alignment Filtering - Usage Guide

## Overview
Filter alignments by mapping quality, flags, regions, and other criteria to prepare clean datasets for downstream analysis.

## Prerequisites
```bash
# samtools
conda install -c bioconda samtools

# pysam
pip install pysam
```

## Quick Start
Tell your AI agent what you want to do:
- "Remove unmapped reads from my BAM file"
- "Filter for high-quality mapped reads (MAPQ >= 30)"
- "Extract reads from specific genomic regions"
- "Remove duplicates and secondary alignments"

## Example Prompts

### Quality Filtering
> "Keep only reads with mapping quality 30 or higher"

> "Remove all unmapped reads from my BAM file"

> "Filter for primary alignments only (no secondary or supplementary)"

### Flag-Based Filtering
> "Remove duplicate reads from my BAM file"

> "Keep only properly paired reads"

> "Extract first-in-pair reads only"

### Region Filtering
> "Extract reads from chr1:1000000-2000000"

> "Get all reads overlapping my target BED file"

> "Extract reads from multiple regions"

### Subsampling
> "Subsample my BAM to 10% of reads"

> "Downsample to approximately 1 million reads"

> "Create a reproducible subset for testing"

## What the Agent Will Do

1. Analyze the current BAM file to understand read composition
2. Apply requested filters using appropriate flags and quality thresholds
3. Count reads before and after filtering
4. Write filtered output to a new BAM file
5. Verify the output file and report filtering statistics

## Understanding SAM FLAGS

| Bit | Flag | Meaning |
|-----|------|---------|
| 0x1 | 1 | Read is paired |
| 0x2 | 2 | Read is properly paired |
| 0x4 | 4 | Read is unmapped |
| 0x8 | 8 | Mate is unmapped |
| 0x10 | 16 | Read is on reverse strand |
| 0x20 | 32 | Mate is on reverse strand |
| 0x40 | 64 | First read in pair |
| 0x80 | 128 | Second read in pair |
| 0x100 | 256 | Secondary alignment |
| 0x200 | 512 | Failed vendor QC |
| 0x400 | 1024 | PCR or optical duplicate |
| 0x800 | 2048 | Supplementary alignment |

### Decoding FLAGS
```bash
samtools flags 99   # 0x63 99 PAIRED,PROPER_PAIR,MREVERSE,READ1
samtools flags 147  # 0x93 147 PAIRED,PROPER_PAIR,REVERSE,READ2
```

## Common Filter Patterns

### -f FLAG: Require These Bits
```bash
samtools view -f 1 input.bam    # Only paired reads
samtools view -f 2 input.bam    # Only properly paired
samtools view -f 3 input.bam    # Paired AND properly paired
```

### -F FLAG: Exclude These Bits
```bash
samtools view -F 4 input.bam      # Remove unmapped
samtools view -F 1024 input.bam   # Remove duplicates
samtools view -F 2304 input.bam   # Remove secondary and supplementary
```

### Combined Filters
```bash
# Properly paired AND not duplicates
samtools view -f 2 -F 1024 input.bam

# Primary mapped reads (most common)
samtools view -F 2308 input.bam  # 2308 = 4 + 256 + 2048

# Variant calling prep
samtools view -f 2 -F 3332 -q 20 -o clean.bam input.bam
```

## Common Commands

### Basic Quality Filter
```bash
samtools view -F 4 -o mapped.bam input.bam              # Keep mapped only
samtools view -F 1024 -o nodup.bam input.bam            # Remove duplicates
samtools view -F 2308 -o primary.bam input.bam          # Primary alignments only
samtools view -F 2308 -q 30 -o filtered.bam input.bam   # Standard quality filter
```

### Region Filtering
```bash
samtools view input.bam chr1:1000000-2000000 -o region.bam
samtools view input.bam chr1:1000-2000 chr2:3000-4000 -o regions.bam
samtools view -L targets.bed input.bam -o targets.bam
samtools view -L targets.bed -F 2308 -q 30 -o filtered_targets.bam input.bam
```

### Subsampling
```bash
samtools view -s 0.1 -o subset.bam input.bam      # ~10% of reads
samtools view -s 42.1 -o subset.bam input.bam     # Reproducible (seed 42)

# Downsample to target count
total=$(samtools view -c input.bam)
frac=$(echo "scale=6; 1000000 / $total" | bc)
samtools view -s "$frac" -o subset.bam input.bam
```

### Output Options
```bash
samtools view -b -F 4 -o output.bam input.bam              # Output BAM
samtools view -C -T reference.fa -F 4 -o output.cram input.bam  # Output CRAM
samtools view -c -F 4 input.bam                            # Count only
samtools view -h -F 4 input.bam > output.sam               # With header
```

## Python with pysam

### Simple Filter
```python
import pysam

with pysam.AlignmentFile('input.bam', 'rb') as infile:
    with pysam.AlignmentFile('filtered.bam', 'wb', header=infile.header) as outfile:
        for read in infile:
            if read.is_unmapped:
                continue
            if read.mapping_quality < 30:
                continue
            outfile.write(read)
```

### Configurable Filter
```python
import pysam

class AlignmentFilter:
    def __init__(self, min_mapq=0, remove_duplicates=True, primary_only=True):
        self.min_mapq = min_mapq
        self.remove_duplicates = remove_duplicates
        self.primary_only = primary_only

    def passes(self, read):
        if read.is_unmapped:
            return False
        if read.mapping_quality < self.min_mapq:
            return False
        if self.remove_duplicates and read.is_duplicate:
            return False
        if self.primary_only and (read.is_secondary or read.is_supplementary):
            return False
        return True

filt = AlignmentFilter(min_mapq=30)
with pysam.AlignmentFile('input.bam', 'rb') as infile:
    with pysam.AlignmentFile('filtered.bam', 'wb', header=infile.header) as outfile:
        for read in infile:
            if filt.passes(read):
                outfile.write(read)
```

### Count with Filter
```python
import pysam

def count_with_filter(bam_path, mapq_min=0, exclude_flags=0):
    count = 0
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for read in bam:
            if read.flag & exclude_flags:
                continue
            if read.mapping_quality >= mapq_min:
                count += 1
    return count

before = count_with_filter('input.bam')
after = count_with_filter('input.bam', mapq_min=30, exclude_flags=2308)
print(f'Before: {before}, After: {after}, Removed: {before - after}')
```

## Mapping Quality (MAPQ)

| MAPQ | Error Probability | Meaning |
|------|-------------------|---------|
| 0 | >10% | Multi-mapped or very uncertain |
| 10 | 10% | Low confidence |
| 20 | 1% | Moderate confidence |
| 30 | 0.1% | High confidence |
| 40 | 0.01% | Very high confidence |
| 60 | ~0% | Maximum (BWA) |

## Troubleshooting

### Index Required for Region Filtering
```bash
samtools index input.bam
samtools view input.bam chr1:1000-2000
```

### Check Filter Effect First
```bash
samtools view -c input.bam                      # Before
samtools view -c -F 2308 -q 30 input.bam        # After
samtools view -F 2308 -q 30 -o filtered.bam input.bam  # Apply if reasonable
```

### Verify Output
```bash
samtools quickcheck filtered.bam && echo "OK" || echo "CORRUPT"
samtools flagstat filtered.bam
```

## Tips
- Always count reads before and after filtering to verify expected behavior
- Use `-F 2308` (unmapped + secondary + supplementary) for most downstream analyses
- MAPQ 30 is a good default threshold for high-confidence alignments
- Region filtering requires an index - index first if needed
- Use `-s SEED.FRACTION` format for reproducible subsampling
- Combine `-f` and `-F` flags to build complex filters
- Save filtering commands for reproducibility

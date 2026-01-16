# Alignment Filtering Usage Guide

This guide covers filtering alignments by various criteria.

## Prerequisites

- samtools installed (`conda install -c bioconda samtools`)
- pysam installed (`pip install pysam`)

## Understanding SAM FLAGS

Every alignment has a FLAG field encoding its properties as a bitwise combination:

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
# What does FLAG 99 mean?
samtools flags 99
# 0x63  99  PAIRED,PROPER_PAIR,MREVERSE,READ1

# What does FLAG 147 mean?
samtools flags 147
# 0x93  147  PAIRED,PROPER_PAIR,REVERSE,READ2
```

### Common FLAG Combinations

| FLAG | Hex | Meaning |
|------|-----|---------|
| 99 | 0x63 | Read1, proper pair, mate reverse |
| 147 | 0x93 | Read2, proper pair, on reverse strand |
| 83 | 0x53 | Read1, proper pair, on reverse strand |
| 163 | 0xA3 | Read2, proper pair, mate reverse |
| 4 | 0x4 | Unmapped |
| 256 | 0x100 | Secondary alignment |
| 2048 | 0x800 | Supplementary (chimeric) |

## The -f and -F Options

### -f FLAG: Require bits

Keep reads where ALL specified bits are set:

```bash
# Keep only paired reads (bit 0x1 must be set)
samtools view -f 1 input.bam

# Keep only properly paired (bits 0x1 and 0x2)
samtools view -f 3 input.bam
```

### -F FLAG: Exclude bits

Remove reads where ANY specified bits are set:

```bash
# Remove unmapped (exclude if bit 0x4 is set)
samtools view -F 4 input.bam

# Remove secondary and supplementary (exclude if 0x100 or 0x800)
samtools view -F 2304 input.bam
```

### Combining -f and -F

```bash
# Properly paired AND not duplicates
samtools view -f 2 -F 1024 input.bam

# Primary reads that are mapped
samtools view -F 2308 input.bam
# 2308 = 4 + 256 + 2048 (unmapped, secondary, supplementary)
```

## Common Filter Scenarios

### Keep Only Mapped Reads

```bash
samtools view -F 4 -o mapped.bam input.bam
```

Removes reads with FLAG bit 4 (unmapped).

### Remove Duplicates

```bash
samtools view -F 1024 -o nodup.bam input.bam
```

After running `samtools markdup`, removes marked duplicates.

### Primary Alignments Only

```bash
samtools view -F 2304 -o primary.bam input.bam
```

Removes secondary (256) and supplementary (2048) alignments.

### Standard Quality Filter

```bash
samtools view -F 2308 -q 30 -o filtered.bam input.bam
```

- `-F 2308`: Remove unmapped, secondary, supplementary
- `-q 30`: Require MAPQ >= 30

### Variant Calling Prep

```bash
samtools view -f 2 -F 3332 -q 20 -o clean.bam input.bam
```

- `-f 2`: Require properly paired
- `-F 3332`: Remove unmapped (4), secondary (256), duplicate (1024), supplementary (2048)
- `-q 20`: MAPQ >= 20

## Filtering by Mapping Quality

### What is MAPQ?

MAPQ (Mapping Quality) indicates confidence that the read is mapped correctly:

```
MAPQ = -10 * log10(P(mapping position is wrong))
```

| MAPQ | Error Probability | Meaning |
|------|-------------------|---------|
| 0 | >10% | Multi-mapped or very uncertain |
| 10 | 10% | Low confidence |
| 20 | 1% | Moderate confidence |
| 30 | 0.1% | High confidence |
| 40 | 0.01% | Very high confidence |
| 60 | ~0% | Maximum (BWA) |

### Filter by MAPQ

```bash
# Minimum MAPQ 30
samtools view -q 30 -o highqual.bam input.bam

# Combined with other filters
samtools view -F 4 -q 30 -o filtered.bam input.bam
```

## Filtering by Region

### Single Region

```bash
# Requires index
samtools view input.bam chr1:1000000-2000000 -o region.bam
```

### Multiple Regions

```bash
samtools view input.bam chr1:1000-2000 chr2:3000-4000 chr3:5000-6000 -o regions.bam
```

### BED File Regions

```bash
# targets.bed format: chr start end
samtools view -L targets.bed input.bam -o targets.bam
```

### Combine with Quality Filter

```bash
samtools view -L targets.bed -F 2308 -q 30 -o filtered_targets.bam input.bam
```

## Subsampling

### Random Fraction

```bash
# Keep ~10% of reads
samtools view -s 0.1 -o subset.bam input.bam
```

### Reproducible Subsampling

```bash
# Seed ensures same reads each time
samtools view -s 42.1 -o subset.bam input.bam
# The integer part (42) is the seed
```

### Downsample to Target

```bash
# Downsample to ~1 million reads
total=$(samtools view -c input.bam)
frac=$(echo "scale=6; 1000000 / $total" | bc)
samtools view -s "$frac" -o subset.bam input.bam
```

## Output Options

### Output BAM (Default)

```bash
samtools view -b -F 4 -o output.bam input.bam
```

### Output CRAM

```bash
samtools view -C -T reference.fa -F 4 -o output.cram input.bam
```

### Count Only

```bash
samtools view -c -F 4 input.bam
```

### With Header

```bash
samtools view -h -F 4 input.bam > output.sam
```

## Working with pysam

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

### Count Filtered Reads

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

total = count_with_filter('input.bam')
filtered = count_with_filter('input.bam', mapq_min=30, exclude_flags=2308)
print(f'Before: {total}, After: {filtered}, Removed: {total - filtered}')
```

## Troubleshooting

### Index Required for Region Filtering

```bash
# Error: [E::idx_find_and_load] could not retrieve index
samtools index input.bam
samtools view input.bam chr1:1000-2000
```

### Check Filter Effect Before Applying

```bash
# Count before filtering
samtools view -c input.bam
# 10000000

# Count after filtering
samtools view -c -F 2308 -q 30 input.bam
# 8500000

# Apply if reasonable
samtools view -F 2308 -q 30 -o filtered.bam input.bam
```

### Verify BAM After Filtering

```bash
samtools quickcheck filtered.bam && echo "OK" || echo "CORRUPT"
samtools flagstat filtered.bam
```

## See Also

- [samtools view documentation](http://www.htslib.org/doc/samtools-view.html)
- [SAM FLAG explanation](https://broadinstitute.github.io/picard/explain-flags.html)

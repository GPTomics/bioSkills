# Alignment Sorting - Usage Guide

## Overview
Sort BAM files by coordinate or read name for downstream analysis tools that require specific sort orders.

## Prerequisites
```bash
# samtools
conda install -c bioconda samtools

# pysam
pip install pysam
```

## Quick Start
Tell your AI agent what you want to do:
- "Sort my BAM file by coordinate"
- "Sort reads by name for duplicate marking"
- "Check the current sort order of my BAM file"
- "Sort directly from aligner output in a pipeline"

## Example Prompts

### Coordinate Sorting
> "Sort my BAM file by genomic position"

> "Sort the output from BWA alignment"

> "Verify my BAM is coordinate-sorted"

### Name Sorting
> "Sort my BAM file by read name"

> "Prepare BAM for the fixmate step"

> "Sort for extracting paired FASTQ files"

### Pipeline Operations
> "Align with BWA and sort in one command"

> "Run the complete duplicate marking workflow"

> "Re-sort an incorrectly sorted BAM file"

### Performance Optimization
> "Sort with 8 threads and 4GB memory per thread"

> "Use a fast SSD for temporary files during sorting"

## What the Agent Will Do

1. Check the current sort order of the input BAM
2. Select appropriate sort method (coordinate or name)
3. Configure memory and thread settings for optimal performance
4. Execute the sort operation
5. Verify the output file was created successfully
6. Index the sorted file if coordinate-sorted

## Sort Order Requirements

| Sort Order | Required For |
|------------|--------------|
| Coordinate | Indexing, IGV, variant calling, mpileup |
| Name | fixmate, markdup (first pass), paired FASTQ extraction |

## Common Commands

### Coordinate Sort
```bash
samtools sort -o sorted.bam unsorted.bam

# From aligner pipeline
bwa mem ref.fa reads.fq | samtools sort -o aligned.bam
```

### Name Sort
```bash
samtools sort -n -o namesorted.bam input.bam
```

### Verify Sort Order
```bash
samtools view -H sorted.bam | grep "^@HD"
# @HD VN:1.6 SO:coordinate
```

### Performance Options
```bash
samtools sort -@ 8 -o sorted.bam input.bam        # 8 threads
samtools sort -@ 4 -m 4G -o sorted.bam input.bam  # 4GB per thread
samtools sort -T /scratch/tmp -o sorted.bam input.bam  # Fast temp disk
samtools sort -l 1 -o sorted.bam input.bam        # Fast compression
```

### Duplicate Marking Workflow
```bash
samtools sort -n -@ 4 input.bam | \
    samtools fixmate -m -@ 4 - - | \
    samtools sort -@ 4 - | \
    samtools markdup -@ 4 - final.bam
samtools index final.bam
```

### Collate vs Sort
```bash
# Collate: fast grouping without full sort
samtools collate -u -O input.bam /tmp/prefix | \
    samtools fastq -1 R1.fq -2 R2.fq -0 /dev/null -s /dev/null -

# Full name sort when strict ordering is needed
samtools sort -n -o namesorted.bam input.bam
```

## Python with pysam

### Basic Sort
```python
import pysam

pysam.sort('-o', 'sorted.bam', 'input.bam')
pysam.sort('-n', '-o', 'namesorted.bam', 'input.bam')
pysam.sort('-@', '4', '-m', '2G', '-o', 'sorted.bam', 'input.bam')
```

### Check Sort Order
```python
import pysam

def get_sort_order(bam_path):
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        hd = bam.header.get('HD', {})
        return hd.get('SO', 'unknown')

order = get_sort_order('input.bam')
print(f'Sort order: {order}')
```

### Sort if Needed
```python
import pysam

def ensure_coordinate_sorted(input_bam, output_bam):
    order = get_sort_order(input_bam)
    if order == 'coordinate':
        return input_bam
    pysam.sort('-o', output_bam, input_bam)
    return output_bam
```

## Troubleshooting

### Out of Memory
Reduce per-thread memory or number of threads:
```bash
samtools sort -@ 4 -m 2G -o sorted.bam input.bam
```

### Disk Full
Temp files filling disk. Use different location:
```bash
samtools sort -T /other/disk/tmp -o sorted.bam input.bam
```

### Slow Sorting
- Increase threads: `-@ 8`
- Reduce compression: `-l 1`
- Use fast disk for temp: `-T /ssd/tmp`
- Increase memory: `-m 4G`

### Corrupted Output
If sort is interrupted, output may be truncated. Re-run from original input.

## Tips
- The `-@` flag specifies additional threads (so `-@ 8` uses 9 total)
- Total memory is approximately threads x memory-per-thread
- Coordinate sort then index is the most common workflow
- Name sort is only needed for specific tools like fixmate
- Collate is faster than name sort when you just need pairs together
- Lower compression (`-l 1`) speeds up sorting but increases file size
- Keep original unsorted file until you verify the sorted output

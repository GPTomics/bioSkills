# Alignment Sorting Usage Guide

This guide covers sorting BAM files by coordinate or read name.

## Prerequisites

- samtools installed (`conda install -c bioconda samtools`)
- pysam installed (`pip install pysam`)

## Why Sort?

Different downstream tools require different sort orders:

| Sort Order | Required For |
|------------|--------------|
| Coordinate | Indexing, IGV, variant calling, mpileup |
| Name | fixmate, markdup (first pass), paired FASTQ extraction |

## Coordinate Sorting

Sorts reads by chromosome and position. This is the most common sort order.

### Basic Sort

```bash
samtools sort -o sorted.bam unsorted.bam
```

### From Aligner Output

Most aligners output unsorted SAM. Sort in a pipeline:

```bash
# BWA + sort
bwa mem ref.fa reads.fq | samtools sort -o aligned.bam

# Bowtie2 + sort
bowtie2 -x index -U reads.fq | samtools sort -o aligned.bam

# STAR outputs sorted by default with --outSAMtype BAM SortedByCoordinate
```

### Verify Sort Order

```bash
samtools view -H sorted.bam | grep "^@HD"
# @HD VN:1.6 SO:coordinate
```

## Name Sorting

Sorts reads by query name, grouping paired reads together.

```bash
samtools sort -n -o namesorted.bam input.bam
```

### When to Use Name Sort

1. **Before fixmate**: fixmate requires name-sorted input
2. **Before markdup**: The name-sort -> fixmate -> coord-sort -> markdup workflow
3. **Extracting paired FASTQ**: Ensures R1 and R2 are adjacent

## The Duplicate Marking Workflow

The standard workflow for marking PCR duplicates:

```bash
# 1. Sort by name
samtools sort -n -o step1_namesort.bam input.bam

# 2. Add mate information
samtools fixmate -m step1_namesort.bam step2_fixmate.bam

# 3. Sort by coordinate
samtools sort -o step3_coordsort.bam step2_fixmate.bam

# 4. Mark duplicates
samtools markdup step3_coordsort.bam final.bam

# 5. Index
samtools index final.bam
```

Or as a pipeline:

```bash
samtools sort -n input.bam | \
    samtools fixmate -m - - | \
    samtools sort - | \
    samtools markdup - final.bam

samtools index final.bam
```

## Collate vs Sort

`samtools collate` groups paired reads without full sorting. It's faster when you just need pairs together.

### Use collate for FASTQ extraction

```bash
# Fast paired FASTQ extraction
samtools collate -u -O input.bam /tmp/prefix | \
    samtools fastq -1 R1.fq -2 R2.fq -0 /dev/null -s /dev/null -
```

### Use sort -n when order matters

```bash
# Full name sort when downstream tool needs strict ordering
samtools sort -n -o namesorted.bam input.bam
```

## Performance Optimization

### Multi-threading

Sort is I/O and CPU intensive. Use multiple threads:

```bash
samtools sort -@ 8 -o sorted.bam input.bam
```

The `-@` flag specifies *additional* threads (so `-@ 8` uses 9 total).

### Memory Control

Control memory per thread to avoid OOM errors:

```bash
# 4GB per thread
samtools sort -@ 4 -m 4G -o sorted.bam input.bam
```

Total memory ≈ threads × memory-per-thread.

### Temporary Files

Sort uses temp files for large datasets. Control location:

```bash
# Use fast SSD for temp files
samtools sort -T /scratch/sort_temp -o sorted.bam input.bam
```

### Compression Level

Lower compression = faster but larger files:

```bash
# Fast sorting, higher disk usage
samtools sort -l 1 -o sorted.bam input.bam

# Maximum compression, slower
samtools sort -l 9 -o sorted.bam input.bam
```

## Working with pysam

### Basic Sort

```python
import pysam

pysam.sort('-o', 'sorted.bam', 'input.bam')
```

### Sort by Name

```python
pysam.sort('-n', '-o', 'namesorted.bam', 'input.bam')
```

### Sort with Options

```python
pysam.sort('-@', '4', '-m', '2G', '-o', 'sorted.bam', 'input.bam')
```

### Check Current Sort Order

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

## Common Workflows

### Align, Sort, Index Pipeline

```bash
#!/bin/bash
REF=$1
R1=$2
R2=$3
OUT=$4

bwa mem -t 8 "$REF" "$R1" "$R2" | \
    samtools sort -@ 4 -o "$OUT"

samtools index "$OUT"
```

### Re-sort Existing BAM

```bash
# Check current order
samtools view -H input.bam | grep "^@HD"

# Re-sort if needed
samtools sort -o resorted.bam input.bam
samtools index resorted.bam
```

### Sort Multiple Files

```bash
for bam in *.bam; do
    out="sorted_${bam}"
    samtools sort -@ 4 -o "$out" "$bam"
    samtools index "$out"
done
```

## Troubleshooting

### Out of Memory

Reduce per-thread memory or number of threads:

```bash
# Instead of 8 threads with default memory
samtools sort -@ 4 -m 2G -o sorted.bam input.bam
```

### Disk Full

Temp files are filling disk. Use different location:

```bash
samtools sort -T /other/disk/tmp -o sorted.bam input.bam
```

### Slow Sorting

1. Increase threads: `-@ 8`
2. Reduce compression: `-l 1`
3. Use fast disk for temp: `-T /ssd/tmp`
4. Increase memory: `-m 4G`

### Corrupted Output

If sort is interrupted, output may be truncated. Always re-run from original input.

## See Also

- [samtools sort documentation](http://www.htslib.org/doc/samtools-sort.html)
- [samtools collate documentation](http://www.htslib.org/doc/samtools-collate.html)

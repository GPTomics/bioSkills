# Duplicate Handling Usage Guide

This guide covers marking and removing PCR/optical duplicates from alignment files.

## Prerequisites

- samtools installed (`conda install -c bioconda samtools`)
- pysam installed (`pip install pysam`)

## What Are Duplicates?

### PCR Duplicates

During library preparation, PCR amplification creates multiple copies of some DNA fragments. When sequenced, these produce identical reads that don't represent independent observations.

### Optical Duplicates

On Illumina flowcells, clusters can sometimes be read more than once if they're very close together, creating artificial duplicates.

### Why Mark Them?

- **Variant calling**: Duplicates inflate coverage and can bias allele frequencies
- **ChIP-seq**: Duplicates can create false peaks
- **RNA-seq**: Usually NOT removed (biological duplicates expected for highly-expressed genes)

## The Standard Workflow

samtools duplicate marking requires a specific order of operations:

```
Input BAM (any order)
    │
    ▼
samtools sort -n (name sort)
    │
    ▼
samtools fixmate -m (add mate info)
    │
    ▼
samtools sort (coordinate sort)
    │
    ▼
samtools markdup (mark/remove duplicates)
    │
    ▼
samtools index (create index)
    │
    ▼
Final BAM (coordinate sorted, duplicates marked)
```

### Why This Order?

1. **fixmate** needs reads grouped by name to find mate pairs
2. **markdup** needs coordinate-sorted reads to find duplicates
3. The MC (mate CIGAR) and ms (mate score) tags added by fixmate are used by markdup

## Step-by-Step Commands

### Individual Steps

```bash
# Step 1: Sort by read name
samtools sort -n -@ 4 -o step1_namesort.bam input.bam

# Step 2: Add mate information
samtools fixmate -m -@ 4 step1_namesort.bam step2_fixmate.bam

# Step 3: Sort by coordinate
samtools sort -@ 4 -o step3_coordsort.bam step2_fixmate.bam

# Step 4: Mark duplicates
samtools markdup -@ 4 step3_coordsort.bam step4_marked.bam

# Step 5: Index
samtools index step4_marked.bam

# Clean up intermediate files
rm step1_namesort.bam step2_fixmate.bam step3_coordsort.bam
```

### Pipeline (No Intermediate Files)

```bash
samtools sort -n -@ 4 input.bam | \
    samtools fixmate -m -@ 4 - - | \
    samtools sort -@ 4 - | \
    samtools markdup -@ 4 - marked.bam

samtools index marked.bam
```

### Full Pipeline from Alignment

```bash
bwa mem -t 8 reference.fa R1.fq R2.fq | \
    samtools sort -n -@ 4 - | \
    samtools fixmate -m -@ 4 - - | \
    samtools sort -@ 4 - | \
    samtools markdup -@ 4 - final.bam

samtools index final.bam
```

## markdup Options

### Mark vs Remove

```bash
# Mark duplicates (FLAG 0x400 set, reads kept)
samtools markdup input.bam marked.bam

# Remove duplicates (duplicates not written to output)
samtools markdup -r input.bam deduped.bam
```

### Statistics Output

```bash
# Print stats to stderr
samtools markdup -s input.bam marked.bam 2> stats.txt

# Write stats to file
samtools markdup -f stats.txt input.bam marked.bam
```

Stats file format:
```
COMMAND: samtools markdup ...
READ: 10000000
WRITTEN: 10000000
EXCLUDED: 0
EXAMINED: 9500000
PAIRED: 9500000
SINGLE: 500000
DUPLICATE PAIR: 500000
DUPLICATE SINGLE: 25000
DUPLICATE PAIR OPTICAL: 10000
DUPLICATE SINGLE OPTICAL: 500
DUPLICATE TOTAL: 525000
ESTIMATED_LIBRARY_SIZE: 45000000
```

### Optical Duplicate Detection

```bash
# Increase pixel distance for patterned flowcells (NovaSeq, NextSeq)
samtools markdup -d 2500 input.bam marked.bam

# Default is 100 pixels (appropriate for older sequencers)
```

## Checking Duplicate Rates

### Using flagstat

```bash
samtools flagstat marked.bam
```

Output includes:
```
10000000 + 0 in total
9800000 + 0 primary
200000 + 0 secondary
0 + 0 supplementary
525000 + 0 duplicates        # <-- This line
```

### Calculate Percentage

```bash
total=$(samtools view -c -F 256 marked.bam)  # Primary alignments
dups=$(samtools view -c -f 1024 -F 256 marked.bam)  # Primary duplicates
echo "Duplicate rate: $(echo "scale=2; $dups * 100 / $total" | bc)%"
```

### Quick Duplicate Count

```bash
# Total duplicates
samtools view -c -f 1024 marked.bam

# Non-duplicates
samtools view -c -F 1024 marked.bam
```

## Working with pysam

### Full Workflow

```python
import pysam
import os

def mark_duplicates(input_bam, output_bam, threads=4):
    temp_ns = 'temp_namesort.bam'
    temp_fm = 'temp_fixmate.bam'
    temp_cs = 'temp_coordsort.bam'

    try:
        pysam.sort('-n', '-@', str(threads), '-o', temp_ns, input_bam)
        pysam.fixmate('-m', '-@', str(threads), temp_ns, temp_fm)
        pysam.sort('-@', str(threads), '-o', temp_cs, temp_fm)
        pysam.markdup('-@', str(threads), temp_cs, output_bam)
        pysam.index(output_bam)
    finally:
        for f in [temp_ns, temp_fm, temp_cs]:
            if os.path.exists(f):
                os.remove(f)

mark_duplicates('input.bam', 'marked.bam')
```

### Check Duplicate Rate

```python
import pysam

def duplicate_rate(bam_path):
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        total = 0
        duplicates = 0
        for read in bam:
            if read.is_secondary or read.is_supplementary:
                continue
            total += 1
            if read.is_duplicate:
                duplicates += 1

    rate = duplicates / total * 100 if total > 0 else 0
    return {'total': total, 'duplicates': duplicates, 'rate': rate}

stats = duplicate_rate('marked.bam')
print(f'Total: {stats["total"]:,}')
print(f'Duplicates: {stats["duplicates"]:,}')
print(f'Rate: {stats["rate"]:.2f}%')
```

### Filter Duplicates

```python
import pysam

def remove_duplicates(input_bam, output_bam):
    with pysam.AlignmentFile(input_bam, 'rb') as infile:
        with pysam.AlignmentFile(output_bam, 'wb', header=infile.header) as outfile:
            for read in infile:
                if not read.is_duplicate:
                    outfile.write(read)

remove_duplicates('marked.bam', 'nodup.bam')
```

## Expected Duplicate Rates

| Library Type | Typical Rate | Concerning |
|--------------|--------------|------------|
| WGS | 5-15% | >25% |
| Exome | 10-20% | >40% |
| ChIP-seq | 10-30% | >50% |
| ATAC-seq | 20-40% | >60% |
| Targeted panels | 30-60% | >80% |

High duplicate rates suggest:
- Low library complexity
- Over-amplification during PCR
- Low input DNA quantity
- Need to sequence more deeply for unique coverage

## Troubleshooting

### "mate not found"

Input to fixmate is not name-sorted:

```bash
samtools sort -n -o namesorted.bam input.bam
samtools fixmate -m namesorted.bam fixmate.bam
```

### "no MC tag found"

fixmate was not run with `-m` flag:

```bash
samtools fixmate -m namesorted.bam fixmate.bam  # Include -m
```

### "not coordinate sorted"

Input to markdup is not coordinate-sorted:

```bash
samtools sort -o coordsorted.bam fixmate.bam
samtools markdup coordsorted.bam marked.bam
```

### High Memory Usage

markdup loads read information into memory. For very large files:

```bash
# Increase threads to parallelize
samtools markdup -@ 8 input.bam marked.bam
```

## Alternatives

### samblaster (Inline with Alignment)

```bash
bwa mem ref.fa R1.fq R2.fq | samblaster | samtools sort -o marked.bam
```

### Picard MarkDuplicates

```bash
java -jar picard.jar MarkDuplicates \
    I=input.bam \
    O=marked.bam \
    M=metrics.txt \
    REMOVE_DUPLICATES=false
```

## See Also

- [samtools markdup documentation](http://www.htslib.org/doc/samtools-markdup.html)
- [samtools fixmate documentation](http://www.htslib.org/doc/samtools-fixmate.html)

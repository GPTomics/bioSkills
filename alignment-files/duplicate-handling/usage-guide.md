# Duplicate Handling - Usage Guide

## Overview
Mark and remove PCR and optical duplicates from alignment files to prevent bias in downstream variant calling and peak detection.

## Prerequisites
```bash
# samtools
conda install -c bioconda samtools

# pysam
pip install pysam
```

## Quick Start
Tell your AI agent what you want to do:
- "Mark duplicates in my BAM file"
- "Remove PCR duplicates from my alignment"
- "Calculate the duplicate rate for my sample"
- "Run the full duplicate marking workflow"

## Example Prompts

### Marking Duplicates
> "Mark duplicates in sample.bam and keep them flagged"

> "Remove duplicates completely from my BAM file"

> "Mark duplicates and generate statistics"

### Workflow Operations
> "Run the complete fixmate and markdup pipeline"

> "Prepare my BAM for duplicate marking with fixmate"

> "Mark duplicates with optical duplicate detection for NovaSeq data"

### Quality Assessment
> "Calculate the duplicate rate for my sample"

> "Check if my duplicate rate is acceptable for WGS"

> "Compare duplicate rates across my samples"

## What the Agent Will Do

1. Sort the BAM file by read name if not already sorted
2. Add mate information using fixmate with the -m flag
3. Re-sort by coordinate position
4. Mark (or remove) duplicates using markdup
5. Generate duplicate statistics
6. Index the final BAM file
7. Report the duplicate rate and assess quality

## The Standard Workflow

```
Input BAM (any order)
    |
    v
samtools sort -n (name sort)
    |
    v
samtools fixmate -m (add mate info)
    |
    v
samtools sort (coordinate sort)
    |
    v
samtools markdup (mark/remove duplicates)
    |
    v
samtools index (create index)
    |
    v
Final BAM (coordinate sorted, duplicates marked)
```

## Common Commands

### Step-by-Step
```bash
samtools sort -n -@ 4 -o namesort.bam input.bam
samtools fixmate -m -@ 4 namesort.bam fixmate.bam
samtools sort -@ 4 -o coordsort.bam fixmate.bam
samtools markdup -@ 4 coordsort.bam marked.bam
samtools index marked.bam
```

### Pipeline (No Intermediate Files)
```bash
samtools sort -n -@ 4 input.bam | \
    samtools fixmate -m -@ 4 - - | \
    samtools sort -@ 4 - | \
    samtools markdup -@ 4 - marked.bam
samtools index marked.bam
```

### Mark vs Remove
```bash
samtools markdup input.bam marked.bam           # Mark only (FLAG 0x400)
samtools markdup -r input.bam deduped.bam       # Remove duplicates
```

### With Statistics
```bash
samtools markdup -s input.bam marked.bam 2> stats.txt
samtools markdup -f stats.txt input.bam marked.bam
```

### Optical Duplicates (Patterned Flowcells)
```bash
samtools markdup -d 2500 input.bam marked.bam  # NovaSeq, NextSeq
```

### Check Duplicate Rate
```bash
samtools flagstat marked.bam | grep duplicates

# Calculate percentage
total=$(samtools view -c -F 256 marked.bam)
dups=$(samtools view -c -f 1024 -F 256 marked.bam)
echo "Duplicate rate: $(echo "scale=2; $dups * 100 / $total" | bc)%"
```

## Python with pysam

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
```

### Check Duplicate Rate
```python
import pysam

def duplicate_rate(bam_path):
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        total, duplicates = 0, 0
        for read in bam:
            if read.is_secondary or read.is_supplementary:
                continue
            total += 1
            if read.is_duplicate:
                duplicates += 1
    return {'total': total, 'duplicates': duplicates, 'rate': duplicates / total * 100 if total > 0 else 0}
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
markdup loads read information into memory. Increase parallelization:
```bash
samtools markdup -@ 8 input.bam marked.bam
```

## Tips
- Always use the full workflow: name-sort, fixmate -m, coord-sort, markdup
- Mark duplicates before variant calling to prevent bias
- For RNA-seq, usually do NOT remove duplicates (biological duplicates expected)
- The -m flag in fixmate is essential for markdup to work properly
- Use -d 2500 for patterned flowcells (NovaSeq, NextSeq)
- Keep the marked BAM even if duplicate rate seems high - filtering can be done later
- Use samblaster for inline marking during alignment

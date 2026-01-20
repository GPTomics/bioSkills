# Pileup Generation - Usage Guide

## Overview
Generate pileup data showing all reads covering each genomic position for variant calling and position-level analysis.

## Prerequisites
```bash
# samtools and bcftools
conda install -c bioconda samtools bcftools

# pysam
pip install pysam

# Reference must be indexed
samtools faidx reference.fa
```

## Quick Start
Tell your AI agent what you want to do:
- "Generate pileup for variant calling"
- "Count alleles at a specific position"
- "Call variants using bcftools"
- "Check read support at position chr1:1000000"

## Example Prompts

### Basic Pileup
> "Generate text pileup for my BAM file"

> "Create pileup for chr1:1000000-2000000 only"

> "Generate pileup with quality filtering (MAPQ 20, baseQ 20)"

### Variant Calling
> "Call variants using samtools mpileup and bcftools"

> "Generate BCF file for efficient variant calling"

> "Call variants from multiple samples together"

### Position Analysis
> "Count allele frequencies at position chr1:1000000"

> "Find all positions with alternative alleles above 10%"

> "Check read support for known variants"

## What the Agent Will Do

1. Verify the BAM file is indexed and reference is available
2. Apply quality filters (mapping quality, base quality)
3. Generate pileup for the requested regions
4. Parse pileup output to extract allele counts
5. Call variants or report position-level statistics
6. Output results in requested format (text, BCF, VCF)

## Understanding Pileup Format

### Text Pileup Output
```
chr1    1000    A    15    ....,,..,..G...    FFFFFFFGFFFFFFFF
```

| Column | Value | Description |
|--------|-------|-------------|
| 1 | chr1 | Chromosome |
| 2 | 1000 | Position (1-based) |
| 3 | A | Reference base |
| 4 | 15 | Depth |
| 5 | ....,,..,..G... | Read bases |
| 6 | FFF... | Base qualities (ASCII - 33) |

### Read Bases Encoding
| Symbol | Meaning |
|--------|---------|
| `.` | Matches reference (forward strand) |
| `,` | Matches reference (reverse strand) |
| `A/C/G/T` | Mismatch (uppercase = forward) |
| `a/c/g/t` | Mismatch (lowercase = reverse) |
| `*` | Deletion |
| `+NNNbases` | Insertion |
| `-NNNbases` | Deletion |
| `^Q` | Start of read (Q = MAPQ) |
| `$` | End of read |

## Common Commands

### Basic Pileup
```bash
samtools mpileup -f reference.fa input.bam > pileup.txt
```

### Quality Filtering
```bash
samtools mpileup -f reference.fa -q 20 -Q 20 input.bam
```
- `-q 20`: Minimum mapping quality
- `-Q 20`: Minimum base quality

### Region-Specific
```bash
samtools mpileup -f reference.fa -r chr1:1000000-2000000 input.bam
samtools mpileup -f reference.fa -l targets.bed input.bam
```

### Variant Calling Pipeline
```bash
# Simple pipeline
samtools mpileup -f reference.fa input.bam | bcftools call -mv -o variants.vcf

# With quality filtering
samtools mpileup -f reference.fa -q 20 -Q 20 input.bam | \
    bcftools call -mv -Oz -o variants.vcf.gz
bcftools index variants.vcf.gz

# BCF intermediate (more efficient)
samtools mpileup -f reference.fa -g -o raw.bcf input.bam
bcftools call -mv raw.bcf -o variants.vcf

# Multi-sample calling
samtools mpileup -f reference.fa sample1.bam sample2.bam sample3.bam | \
    bcftools call -mv -o variants.vcf
```

### Performance Options
```bash
# Limit max depth (reduces memory)
samtools mpileup -f reference.fa -d 1000 input.bam

# Parallel by chromosome
for chr in chr1 chr2 chr3; do
    samtools mpileup -f reference.fa -r "$chr" input.bam | \
        bcftools call -mv -o "${chr}.vcf" &
done
wait
```

## Python with pysam

### Basic Pileup Iteration
```python
import pysam

with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for pileup_column in bam.pileup('chr1', 1000000, 1001000, truncate=True):
        pos = pileup_column.pos
        depth = pileup_column.n
        print(f'chr1:{pos+1} depth={depth}')
```

### Count Alleles at Position
```python
import pysam
from collections import Counter

def count_alleles(bam_path, chrom, pos, min_qual=20):
    alleles = Counter()
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for pileup_column in bam.pileup(chrom, pos, pos + 1, truncate=True, min_base_quality=min_qual):
            if pileup_column.pos != pos:
                continue
            for pileup_read in pileup_column.pileups:
                if pileup_read.is_del:
                    alleles['DEL'] += 1
                elif pileup_read.is_refskip:
                    pass
                else:
                    qpos = pileup_read.query_position
                    base = pileup_read.alignment.query_sequence[qpos].upper()
                    alleles[base] += 1
    return dict(alleles)

counts = count_alleles('input.bam', 'chr1', 1000000)
total = sum(counts.values())
for base, count in sorted(counts.items(), key=lambda x: -x[1]):
    print(f'{base}: {count} ({count/total*100:.1f}%)')
```

### Find Variants
```python
import pysam
from collections import Counter

def find_variants(bam_path, ref_path, chrom, start, end, min_depth=10, min_alt_freq=0.1):
    variants = []
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        with pysam.FastaFile(ref_path) as ref:
            for pileup_column in bam.pileup(chrom, start, end, truncate=True, min_base_quality=20):
                pos = pileup_column.pos
                ref_base = ref.fetch(chrom, pos, pos + 1).upper()
                alleles = Counter()
                for pileup_read in pileup_column.pileups:
                    if pileup_read.is_del or pileup_read.is_refskip:
                        continue
                    qpos = pileup_read.query_position
                    base = pileup_read.alignment.query_sequence[qpos].upper()
                    alleles[base] += 1
                total = sum(alleles.values())
                if total < min_depth:
                    continue
                for base, count in alleles.items():
                    if base != ref_base:
                        freq = count / total
                        if freq >= min_alt_freq:
                            variants.append({'chrom': chrom, 'pos': pos + 1, 'ref': ref_base, 'alt': base, 'depth': total, 'alt_count': count, 'freq': freq})
    return variants
```

### Access Individual Reads
```python
import pysam

with pysam.AlignmentFile('input.bam', 'rb') as bam:
    for pileup_column in bam.pileup('chr1', 1000000, 1000001, truncate=True):
        print(f'Position {pileup_column.pos}, depth {pileup_column.n}')
        for pileup_read in pileup_column.pileups:
            if pileup_read.is_del:
                print(f'  {pileup_read.alignment.query_name}: DELETION')
                continue
            qpos = pileup_read.query_position
            base = pileup_read.alignment.query_sequence[qpos]
            qual = pileup_read.alignment.query_qualities[qpos]
            strand = '-' if pileup_read.alignment.is_reverse else '+'
            print(f'  {pileup_read.alignment.query_name}: {base} Q{qual} {strand}')
```

## Troubleshooting

### "No sequences in common"
Reference doesn't match the BAM:
```bash
samtools view -H input.bam | grep "^@SQ"  # Check BAM chromosomes
grep "^>" reference.fa                      # Check reference chromosomes
```

### Empty Output
- Check region exists: `samtools view input.bam chr1:1000-2000 | head`
- Check reference is indexed: `ls reference.fa.fai`
- Check chromosome name match (chr1 vs 1)

### Memory Issues
Reduce max depth:
```bash
samtools mpileup -f reference.fa -d 500 input.bam
```

### Slow Processing
- Use BCF output instead of text pileup
- Process chromosomes in parallel
- Filter to target regions with `-l targets.bed`

## Tips
- Always use quality filtering (`-q 20 -Q 20`) for variant calling
- BCF format is much faster to process than text pileup
- Use `truncate=True` in pysam for exact region boundaries
- Multi-sample calling with bcftools is preferred over merging individual VCFs
- Limit max depth with `-d` flag for high-coverage samples
- Text pileup is useful for debugging but inefficient for large-scale analysis
- pysam's pileup uses 0-based coordinates

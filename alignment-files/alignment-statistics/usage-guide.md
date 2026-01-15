# Alignment Statistics Usage Guide

This guide covers generating QC statistics from alignment files.

## Prerequisites

- samtools installed (`conda install -c bioconda samtools`)
- pysam installed (`pip install pysam`)
- For plotting: `conda install -c bioconda gnuplot`

## Choosing the Right Command

| Need | Command | Speed | Notes |
|------|---------|-------|-------|
| Quick read counts | flagstat | <1 sec | Always start here |
| Per-chromosome counts | idxstats | <1 sec | Needs index |
| Full QC report | stats | 1-10 min | Generates plots |
| Coverage summary | coverage | 1-5 min | Per-region stats |
| Per-position depth | depth | 5-30 min | Very detailed |

## flagstat - Quick Overview

The fastest way to get alignment statistics:

```bash
samtools flagstat sample.bam
```

### Understanding the Output

```
10000000 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
50000 + 0 supplementary
500000 + 0 duplicates
9800000 + 0 mapped (98.00% : N/A)
9950000 + 0 paired in sequencing
4975000 + 0 read1
4975000 + 0 read2
9700000 + 0 properly paired (97.49% : N/A)
9750000 + 0 with itself and mate mapped
100000 + 0 singletons (1.01% : N/A)
25000 + 0 with mate mapped to a different chr
10000 + 0 with mate mapped to a different chr (mapQ>=5)
```

Key metrics:
- **mapped**: Percentage of reads that aligned
- **properly paired**: Both mates aligned in expected orientation/distance
- **duplicates**: PCR/optical duplicates (after markdup)
- **singletons**: Paired reads where only one mate mapped

### Parse flagstat

```bash
# Get mapping rate
samtools flagstat sample.bam | grep "mapped" | head -1

# Extract percentage
samtools flagstat sample.bam | awk '/mapped/ && !/mate/ {print $5}'
```

## idxstats - Per-Chromosome Counts

Fast chromosome-level statistics (requires index):

```bash
samtools idxstats sample.bam
```

Output:
```
chr1    248956422    5000000    1000
chr2    242193529    4800000    800
chrM    16569        50000      100
*       0            0          150000
```

Columns: chromosome, length, mapped reads, unmapped reads placed on chromosome

### Common Analyses

```bash
# Total reads per chromosome
samtools idxstats sample.bam | sort -k3 -nr | head

# Mitochondrial fraction (contamination check)
samtools idxstats sample.bam | awk '
    BEGIN {OFS="\t"}
    /^chrM/ {mt = $3}
    {total += $3}
    END {printf "Mitochondrial: %.2f%%\n", mt/total*100}'

# Sex check (X/Y ratio)
samtools idxstats sample.bam | awk '
    /^chrX/ {x = $3}
    /^chrY/ {y = $3}
    END {printf "X:Y ratio = %.2f\n", x/(y+1)}'
```

## stats - Comprehensive Report

Detailed statistics including insert sizes, quality distributions:

```bash
samtools stats sample.bam > sample_stats.txt
```

### Summary Numbers (SN Lines)

```bash
samtools stats sample.bam | grep "^SN" | cut -f 2-
```

Key fields:
- `raw total sequences`: Total reads in file
- `reads mapped`: Number of mapped reads
- `reads properly paired`: Properly paired reads
- `insert size average`: Mean insert size
- `insert size standard deviation`: Insert size spread
- `average length`: Mean read length
- `error rate`: Mismatch rate per base

### Generate Plots

```bash
samtools stats sample.bam > stats.txt
plot-bamstats -p stats_plots/ stats.txt
```

Creates:
- `stats_plots-gc-content.png` - GC content distribution
- `stats_plots-insert-size.png` - Insert size histogram
- `stats_plots-quals.png` - Quality score distribution
- `stats_plots-coverage.png` - Coverage distribution

### Stats for Specific Region

```bash
samtools stats sample.bam chr1:1000000-2000000 > region_stats.txt
```

## coverage - Per-Region Summary

Faster than depth for region-level coverage:

```bash
samtools coverage sample.bam
```

Output:
```
#rname  startpos  endpos    numreads  covbases  coverage  meandepth  meanbaseq  meanmapq
chr1    1         248956422 5000000   240000000 96.40     25.5       35.2       58.1
chr2    1         242193529 4800000   235000000 97.03     24.8       35.1       57.9
```

### Coverage for Target Regions

```bash
# Specific region
samtools coverage -r chr1:1000000-2000000 sample.bam

# From BED file
samtools coverage -b targets.bed sample.bam
```

### Histogram Mode

```bash
samtools coverage -m sample.bam
```

Shows ASCII histogram of coverage distribution.

## depth - Per-Position Detail

Per-base read depth (slow but detailed):

```bash
samtools depth sample.bam > depth.txt
```

### Include Zero-Depth Positions

```bash
samtools depth -a sample.bam > depth_with_zeros.txt
```

### Specific Region

```bash
samtools depth -r chr1:1000000-1001000 sample.bam
```

### Calculate Statistics

```bash
# Mean depth
samtools depth sample.bam | awk '{sum+=$3; n++} END {print "Mean:", sum/n}'

# Median depth (slower, loads into memory)
samtools depth sample.bam | cut -f3 | sort -n | awk '{a[NR]=$1} END {print "Median:", a[int(NR/2)]}'

# Coverage at different thresholds
samtools depth sample.bam | awk '
    {total++; if($3>=10) c10++; if($3>=20) c20++; if($3>=30) c30++}
    END {
        printf ">=10x: %.1f%%\n", c10/total*100
        printf ">=20x: %.1f%%\n", c20/total*100
        printf ">=30x: %.1f%%\n", c30/total*100
    }'
```

### Depth in BED Regions

```bash
samtools depth -b targets.bed sample.bam | awk '
    {sum+=$3; n++}
    END {print "Mean target depth:", sum/n}'
```

## Working with pysam

### Flagstat Equivalent

```python
import pysam

def flagstat(bam_path):
    stats = {'total': 0, 'mapped': 0, 'paired': 0, 'proper': 0, 'duplicate': 0}

    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for read in bam:
            stats['total'] += 1
            if not read.is_unmapped:
                stats['mapped'] += 1
            if read.is_paired:
                stats['paired'] += 1
            if read.is_proper_pair:
                stats['proper'] += 1
            if read.is_duplicate:
                stats['duplicate'] += 1

    stats['map_rate'] = stats['mapped'] / stats['total'] * 100
    stats['proper_rate'] = stats['proper'] / stats['paired'] * 100 if stats['paired'] > 0 else 0

    return stats

stats = flagstat('sample.bam')
for k, v in stats.items():
    print(f'{k}: {v}')
```

### Coverage in Region

```python
import pysam

def region_coverage(bam_path, chrom, start, end):
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        depths = [0] * (end - start)
        for pileup in bam.pileup(chrom, start, end, truncate=True):
            if start <= pileup.pos < end:
                depths[pileup.pos - start] = pileup.n

    covered = sum(1 for d in depths if d > 0)
    mean_depth = sum(depths) / len(depths)

    return {
        'length': len(depths),
        'covered': covered,
        'pct_covered': covered / len(depths) * 100,
        'mean_depth': mean_depth,
        'min_depth': min(depths),
        'max_depth': max(depths)
    }

stats = region_coverage('sample.bam', 'chr1', 1000000, 2000000)
print(f'Coverage: {stats["pct_covered"]:.1f}%')
print(f'Mean depth: {stats["mean_depth"]:.1f}x')
```

## QC Thresholds

### Typical Good Values

| Metric | WGS | Exome | RNA-seq |
|--------|-----|-------|---------|
| Mapping rate | >95% | >95% | >80% |
| Proper pair rate | >90% | >90% | >80% |
| Duplicate rate | <20% | <30% | N/A |
| Mean depth | 30-50x | 50-100x | N/A |
| Target coverage (20x) | >95% | >90% | N/A |

### Red Flags

- Mapping rate <80%: Contamination, wrong reference, poor quality
- High mitochondrial: mtDNA contamination
- Insert size bimodal: Library prep issue
- Error rate >2%: Sequencing quality issue

## Batch Processing

### Process Multiple Files

```bash
for bam in *.bam; do
    echo "=== $bam ===" >> all_stats.txt
    samtools flagstat "$bam" >> all_stats.txt
    echo "" >> all_stats.txt
done
```

### Summary Table

```bash
echo -e "Sample\tTotal\tMapped\tPaired\tDuplicates" > summary.tsv
for bam in *.bam; do
    sample=$(basename "$bam" .bam)
    samtools flagstat "$bam" | awk -v s="$sample" '
        /in total/ {total=$1}
        /mapped \(/ {mapped=$1}
        /properly paired/ {paired=$1}
        /duplicates/ {dup=$1}
        END {print s"\t"total"\t"mapped"\t"paired"\t"dup}
    ' >> summary.tsv
done
```

## See Also

- [samtools stats documentation](http://www.htslib.org/doc/samtools-stats.html)
- [samtools flagstat documentation](http://www.htslib.org/doc/samtools-flagstat.html)
- [samtools depth documentation](http://www.htslib.org/doc/samtools-depth.html)
- [samtools coverage documentation](http://www.htslib.org/doc/samtools-coverage.html)

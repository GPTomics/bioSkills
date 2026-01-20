# Alignment Statistics - Usage Guide

## Overview
Generate QC statistics from alignment files including mapping rates, read counts, coverage depth, and per-chromosome distributions.

## Prerequisites
```bash
# samtools
conda install -c bioconda samtools

# pysam and plotting
pip install pysam matplotlib

# For plot generation
conda install -c bioconda gnuplot
```

## Quick Start
Tell your AI agent what you want to do:
- "Get quick statistics for my BAM file"
- "Calculate coverage depth across my genome"
- "Generate a comprehensive QC report with plots"
- "Check per-chromosome read distribution"

## Example Prompts

### Quick Statistics
> "Run flagstat on my BAM file and interpret the results"

> "How many reads are mapped in sample.bam?"

> "What percentage of reads are properly paired?"

### Per-Chromosome Analysis
> "Show read counts per chromosome using idxstats"

> "Calculate the mitochondrial contamination percentage"

> "Check X/Y ratio for sex determination"

### Coverage Analysis
> "Calculate mean coverage depth for my WGS sample"

> "What percentage of the genome is covered at 20x?"

> "Generate depth statistics for target regions"

### Comprehensive Reports
> "Generate samtools stats and create QC plots"

> "Create a summary table of statistics for all my samples"

## What the Agent Will Do

1. Run appropriate statistics command (flagstat, idxstats, stats, coverage, depth)
2. Parse the output into interpretable metrics
3. Compare metrics against quality thresholds
4. Generate plots if requested
5. Summarize results with pass/fail assessment
6. Provide recommendations for any concerning metrics

## Choosing the Right Command

| Need | Command | Speed |
|------|---------|-------|
| Quick read counts | flagstat | <1 sec |
| Per-chromosome counts | idxstats | <1 sec |
| Full QC report | stats | 1-10 min |
| Coverage summary | coverage | 1-5 min |
| Per-position depth | depth | 5-30 min |

## Common Commands

### flagstat - Quick Overview
```bash
samtools flagstat sample.bam
```

Output interpretation:
- **mapped**: Percentage of reads that aligned
- **properly paired**: Both mates aligned in expected orientation/distance
- **duplicates**: PCR/optical duplicates (after markdup)
- **singletons**: Paired reads where only one mate mapped

### idxstats - Per-Chromosome Counts
```bash
samtools idxstats sample.bam

# Total mapped reads
samtools idxstats sample.bam | awk '{sum += $3} END {print sum}'

# Mitochondrial fraction
samtools idxstats sample.bam | awk '/^chrM/ {mt=$3} {total+=$3} END {printf "MT: %.2f%%\n", mt/total*100}'

# Sex check (X/Y ratio)
samtools idxstats sample.bam | awk '/^chrX/ {x=$3} /^chrY/ {y=$3} END {printf "X:Y = %.2f\n", x/(y+1)}'
```

### stats - Comprehensive Report
```bash
samtools stats sample.bam > sample_stats.txt
grep "^SN" sample_stats.txt | cut -f 2-

# Generate plots
plot-bamstats -p stats_plots/ sample_stats.txt
```

### coverage - Per-Region Summary
```bash
samtools coverage sample.bam
samtools coverage -r chr1:1000000-2000000 sample.bam
samtools coverage -m sample.bam  # ASCII histogram
```

### depth - Per-Position Detail
```bash
samtools depth sample.bam > depth.txt
samtools depth -a sample.bam > depth_with_zeros.txt  # Include zero positions
samtools depth -r chr1:1000000-1001000 sample.bam

# Mean depth
samtools depth sample.bam | awk '{sum+=$3; n++} END {print "Mean:", sum/n}'

# Coverage at thresholds
samtools depth sample.bam | awk '{total++; if($3>=10) c10++; if($3>=20) c20++} END {printf ">=10x: %.1f%% >=20x: %.1f%%\n", c10/total*100, c20/total*100}'
```

## Python with pysam

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
    return stats
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
    return {'length': len(depths), 'covered': covered, 'pct_covered': covered / len(depths) * 100, 'mean_depth': mean_depth}
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

## Troubleshooting

### Missing Index
idxstats requires an index:
```bash
samtools index sample.bam
samtools idxstats sample.bam
```

### Slow depth Command
Use region filtering or sampling:
```bash
samtools depth -r chr1:1-10000000 sample.bam  # Single chromosome
samtools view -s 0.1 sample.bam | samtools depth -  # Sample 10%
```

## Tips
- Start with flagstat for a quick overview before running detailed statistics
- Use idxstats for fast per-chromosome counts without reading alignment data
- samtools stats generates comprehensive metrics but takes longer
- Coverage command is good for region-level summaries
- Depth command provides per-position detail but can be slow for whole genomes
- Always compare metrics across samples to identify outliers
- Keep statistics reports for reproducibility and troubleshooting

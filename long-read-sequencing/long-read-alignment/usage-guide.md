# Long-Read Alignment with minimap2 - Usage Guide

## Overview

minimap2 is the standard aligner for long reads. It's fast, accurate, and handles the higher error rates of long-read data. Different presets optimize for ONT, PacBio, or RNA sequencing.

## When to Use This Skill

- You have Oxford Nanopore FASTQ files
- You have PacBio HiFi or CLR reads
- You need to align long reads for variant calling or assembly

## Installation

```bash
conda install -c bioconda minimap2 samtools
```

## Choosing the Right Preset

| Data Type | Preset |
|-----------|--------|
| Nanopore (any) | map-ont |
| PacBio HiFi (Q20+) | map-hifi |
| PacBio CLR (older) | map-pb |
| cDNA/direct RNA | splice |

## Basic Workflow

```bash
# 1. Align and sort
minimap2 -ax map-ont -t 8 reference.fa reads.fastq.gz | \
    samtools sort -@ 4 -o aligned.bam

# 2. Index BAM
samtools index aligned.bam

# 3. Check alignment rate
samtools flagstat aligned.bam
```

## Best Practices

### Always Use Proper Preset

Wrong preset will give poor alignments. Check your data type.

### Add Read Groups

Required for many downstream tools:
```bash
minimap2 -ax map-ont -R '@RG\tID:run1\tSM:sample1' ref.fa reads.fq.gz
```

### Pre-Index Large References

Saves time for multiple runs:
```bash
minimap2 -d ref.mmi ref.fa
```

## Common Issues

### Low Alignment Rate

- Check correct preset for your data type
- Verify reference matches organism
- Check if reads are very short (< 500bp)

### Slow Performance

- Pre-build index
- Increase threads (-t)
- Use faster storage

## Output

| Extension | Description |
|-----------|-------------|
| .bam | Binary alignment |
| .paf | Pairwise format (no sequence) |

## Resources

- [minimap2 Manual](https://lh3.github.io/minimap2/minimap2.html)
- [minimap2 Paper](https://doi.org/10.1093/bioinformatics/bty191)

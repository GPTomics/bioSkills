# Coverage Analysis - Usage Guide

## Overview

Coverage analysis calculates how well genomic regions are covered by sequencing reads. This skill covers generating bedGraph files with bedtools genomecov and calculating per-feature coverage statistics with bedtools coverage.

## When to Use This Skill

- Generate coverage tracks for genome browsers
- Calculate mean depth across target regions
- Identify regions with low coverage
- Normalize coverage for comparison between samples
- Summarize exon coverage for RNA-seq QC

## Installation

```bash
# bedtools
conda install -c bioconda bedtools

# pybedtools (optional)
pip install pybedtools
```

## Key Concepts

### genomecov vs coverage

| Tool | Purpose | Output |
|------|---------|--------|
| genomecov | Genome-wide coverage | bedGraph or histogram |
| coverage | Coverage per feature | Feature + coverage stats |

### bedGraph Format

bedGraph is a simple 4-column format for continuous data:
```
chr1    0       100     0      # Zero coverage, bases 0-99
chr1    100     200     5.5    # 5.5x coverage, bases 100-199
chr1    200     300     10.2   # 10.2x coverage, bases 200-299
```

## Basic Workflow

### 1. Generate Coverage Track

```bash
# From BAM file
bedtools genomecov -ibam alignments.bam -bg > coverage.bedGraph
```

### 2. Normalize for Comparison

```bash
# Calculate scaling factor (reads per million)
TOTAL=$(samtools view -c alignments.bam)
SCALE=$(echo "scale=10; 1000000/$TOTAL" | bc)

# Generate normalized bedGraph
bedtools genomecov -ibam alignments.bam -bg -scale $SCALE > normalized.bedGraph
```

### 3. Calculate Per-Region Coverage

```bash
# Coverage statistics for target regions
bedtools coverage -a targets.bed -b alignments.bam > target_coverage.bed

# Output columns: original BED + overlaps, covered_bases, region_length, fraction
```

### 4. Find Low-Coverage Regions

```bash
# Regions with <80% coverage
awk '$NF < 0.8' target_coverage.bed > low_coverage.bed
```

## Common Patterns

### QC Target Capture

```bash
# Summarize coverage across targets
bedtools coverage -a targets.bed -b alignments.bam | \
    awk 'BEGIN{sum=0; n=0} {sum+=$NF; n++} END{print "Mean coverage fraction:", sum/n}'
```

### RNA-seq Exon Coverage

```bash
# Handle spliced alignments
bedtools coverage -a exons.bed -b alignments.bam -split > exon_coverage.bed
```

### ChIP-seq Signal

```bash
# Generate signal track for browser
bedtools genomecov -ibam chip.bam -bg > chip_signal.bedGraph
```

## Common Issues

### BAM Not Sorted

**Problem:** "BAM must be sorted"
**Solution:** Sort BAM first

```bash
samtools sort alignments.bam -o sorted.bam
samtools index sorted.bam
bedtools genomecov -ibam sorted.bam -bg > coverage.bedGraph
```

### Missing Genome File

**Problem:** BED input requires -g flag
**Solution:** Provide genome file

```bash
bedtools genomecov -i regions.bed -g genome.txt -bg > coverage.bedGraph
```

### Large Output Files

**Problem:** bedGraph files can be very large
**Solution:** Convert to bigWig for efficiency (see bigwig-tracks skill)

## Resources

- [bedtools genomecov](https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html)
- [bedtools coverage](https://bedtools.readthedocs.io/en/latest/content/tools/coverage.html)
- [bedGraph format](https://genome.ucsc.edu/goldenPath/help/bedgraph.html)

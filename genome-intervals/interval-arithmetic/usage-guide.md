# Interval Arithmetic - Usage Guide

## Overview

Interval arithmetic covers the core set operations on genomic intervals: finding overlaps (intersect), removing overlaps (subtract), combining intervals (merge), and finding uncovered regions (complement). These operations form the foundation of genomic analysis.

## When to Use This Skill

- Find overlapping regions between two sets of intervals
- Remove unwanted regions from your intervals
- Merge overlapping or adjacent intervals
- Find genomic regions not covered by your intervals
- Calculate similarity between interval sets (Jaccard)
- Test statistical significance of overlaps (Fisher)

## Installation

```bash
# bedtools (required)
conda install -c bioconda bedtools

# pybedtools (optional, for Python)
pip install pybedtools
```

## Choosing the Right Operation

| Goal | Operation | Example |
|------|-----------|---------|
| Which A intervals overlap B? | intersect -u | Peaks in promoters |
| Which A intervals don't overlap B? | intersect -v | Peaks outside genes |
| What regions overlap between A and B? | intersect | Overlapping portions |
| Remove B regions from A | subtract | Remove blacklist |
| Combine overlapping intervals | merge | Merge replicates |
| Find uncovered genome | complement | Non-coding regions |

## Basic Workflow

### 1. Prepare Input Files

```bash
# Sort BED files (required for many operations)
bedtools sort -i peaks.bed > peaks.sorted.bed
bedtools sort -i genes.bed > genes.sorted.bed
```

### 2. Find Overlaps

```bash
# Peaks that overlap genes
bedtools intersect -a peaks.sorted.bed -b genes.sorted.bed -u > peaks_in_genes.bed

# Peaks that don't overlap genes
bedtools intersect -a peaks.sorted.bed -b genes.sorted.bed -v > peaks_outside_genes.bed
```

### 3. Merge Replicates

```bash
# Combine peaks from multiple replicates
cat rep1.bed rep2.bed rep3.bed | \
    bedtools sort | \
    bedtools merge -d 100 > consensus_peaks.bed
```

## Common Issues

### Empty Output

**Problem:** No overlaps found
**Solutions:**
- Check chromosome naming (chr1 vs 1)
- Verify coordinate system (0-based vs 1-based)
- Check file format with `head -5`

### Unsorted Input Error

**Problem:** "Input must be sorted"
**Solution:** Sort first

```bash
bedtools sort -i input.bed | bedtools merge > output.bed
```

### Memory Issues

**Problem:** Out of memory on large files
**Solution:** Use streaming or sorted input

```bash
# Stream through operations
bedtools sort -i huge.bed | bedtools merge > merged.bed
```

## Python Integration

```python
import pybedtools
import pandas as pd

# Load and intersect
a = pybedtools.BedTool('peaks.bed')
b = pybedtools.BedTool('genes.bed')
result = a.intersect(b, u=True)

# Convert to DataFrame for analysis
df = result.to_dataframe()
print(f'Found {len(df)} overlapping peaks')

# Cleanup
pybedtools.cleanup()
```

## Performance Tips

1. **Sort once:** Sort input files before multiple operations
2. **Use streams:** Pipe operations together
3. **Index BAM files:** For BAM/BED intersections
4. **Clean up:** Remove temp files with `pybedtools.cleanup()`

## Resources

- [bedtools intersect](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html)
- [bedtools merge](https://bedtools.readthedocs.io/en/latest/content/tools/merge.html)
- [pybedtools tutorial](https://daler.github.io/pybedtools/topical-intersections.html)

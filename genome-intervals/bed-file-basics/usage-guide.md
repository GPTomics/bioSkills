# BED File Basics - Usage Guide

## Overview

BED (Browser Extensible Data) is the standard format for representing genomic intervals. This skill covers creating, reading, validating, and converting BED files using bedtools (CLI) and pybedtools (Python).

## When to Use This Skill

- You need to create BED files from coordinates
- You want to validate BED file format
- You need to convert between BED and other formats
- You want to filter or sort interval files

## Installation

```bash
# bedtools (CLI)
conda install -c bioconda bedtools

# pybedtools (Python)
pip install pybedtools

# Verify
bedtools --version  # Should be 2.31+
python -c "import pybedtools; print(pybedtools.__version__)"
```

## Key Concepts

### Coordinate System

BED uses **0-based, half-open** coordinates:

| System | First base | Interval 100-200 |
|--------|------------|------------------|
| BED (0-based) | 0 | Bases 100-199 |
| GFF/VCF (1-based) | 1 | Bases 100-200 |

When converting:
- BED to GFF: add 1 to start
- GFF to BED: subtract 1 from start

### BED Columns

| Columns | Name | Required Fields |
|---------|------|-----------------|
| BED3 | Minimal | chr, start, end |
| BED4 | Named | + name |
| BED5 | Scored | + score |
| BED6 | Stranded | + strand |
| BED12 | Full | + thick, rgb, blocks |

## Basic Workflow

### 1. Create BED File

```python
import pybedtools
import pandas as pd

# From DataFrame
df = pd.DataFrame({
    'chrom': ['chr1', 'chr1', 'chr2'],
    'start': [100, 500, 200],
    'end': [200, 600, 400],
    'name': ['peak1', 'peak2', 'peak3']
})
bed = pybedtools.BedTool.from_dataframe(df)
bed.saveas('peaks.bed')
```

### 2. Validate

```bash
# Check for issues
awk '$2 >= $3' peaks.bed  # Should be empty
bedtools sort -i peaks.bed > /dev/null  # Should succeed
```

### 3. Sort

```bash
bedtools sort -i peaks.bed > peaks.sorted.bed
```

### 4. Filter

```python
import pybedtools
bed = pybedtools.BedTool('peaks.sorted.bed')
filtered = bed.filter(lambda x: len(x) >= 100)
filtered.saveas('peaks.filtered.bed')
```

## Common Issues

### Wrong Coordinate System

**Problem:** Coordinates off by one
**Solution:** Check if source uses 1-based coordinates (GFF, VCF) and convert

```python
# 1-based to 0-based
df['start'] = df['start'] - 1
```

### Unsorted Input

**Problem:** bedtools operations fail or give wrong results
**Solution:** Sort before operations

```bash
bedtools sort -i input.bed > sorted.bed
```

### Memory Issues with pybedtools

**Problem:** Temp files accumulate
**Solution:** Clean up explicitly

```python
import pybedtools
# At end of script
pybedtools.cleanup()
```

## Resources

- [UCSC BED Format Specification](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)
- [bedtools Documentation](https://bedtools.readthedocs.io/)
- [pybedtools Documentation](https://daler.github.io/pybedtools/)

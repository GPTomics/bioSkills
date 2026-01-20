# Interval Arithmetic - Usage Guide

## Overview
Interval arithmetic covers the core set operations on genomic intervals: finding overlaps (intersect), removing overlaps (subtract), combining intervals (merge), and finding uncovered regions (complement). These operations form the foundation of genomic analysis.

## Prerequisites
```bash
# bedtools (required)
conda install -c bioconda bedtools

# pybedtools (optional, for Python)
pip install pybedtools
```

## Quick Start
Tell your AI agent what you want to do:
- "Find which of my peaks overlap with promoter regions"
- "Remove blacklisted regions from my peak calls"
- "Merge overlapping peaks from multiple replicates"

## Example Prompts

### Finding Overlaps
> "Find peaks from peaks.bed that overlap with promoters.bed"
> "Get the overlapping portions between my ChIP-seq peaks and enhancer annotations"
> "Which of my peaks do NOT overlap with any gene?"

### Removing Regions
> "Remove ENCODE blacklist regions from my peak file"
> "Subtract exons from gene bodies to get introns"

### Merging Intervals
> "Merge overlapping peaks from my three replicates"
> "Combine peaks that are within 100bp of each other"

### Complement and Gaps
> "Find all genomic regions not covered by my annotations"
> "Get the gaps between my exons"

### Statistical Analysis
> "Calculate Jaccard similarity between two peak sets"
> "Test if my peaks significantly overlap with promoters using Fisher's exact test"

## What the Agent Will Do
1. Sort input BED files if needed
2. Execute the appropriate bedtools/pybedtools operation
3. Apply any filtering criteria (distance thresholds, overlap fractions)
4. Save results and report summary statistics

## Choosing the Right Operation

| Goal | Operation | Example |
|------|-----------|---------|
| Which A intervals overlap B? | intersect -u | Peaks in promoters |
| Which A intervals don't overlap B? | intersect -v | Peaks outside genes |
| What regions overlap between A and B? | intersect | Overlapping portions |
| Remove B regions from A | subtract | Remove blacklist |
| Combine overlapping intervals | merge | Merge replicates |
| Find uncovered genome | complement | Non-coding regions |

## Tips
- Sort input files once before running multiple operations
- Use `-u` flag to get unique overlapping intervals (not duplicated per overlap)
- Use `-v` to invert and get non-overlapping intervals
- Pipe operations together for efficiency: `bedtools sort | bedtools merge`
- Use `pybedtools.cleanup()` to remove temp files after Python operations
- Check chromosome naming consistency between files (chr1 vs 1)

## Resources
- [bedtools intersect](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html)
- [bedtools merge](https://bedtools.readthedocs.io/en/latest/content/tools/merge.html)
- [pybedtools tutorial](https://daler.github.io/pybedtools/topical-intersections.html)

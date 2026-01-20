# BED File Basics - Usage Guide

## Overview
BED (Browser Extensible Data) is the standard format for representing genomic intervals. This skill covers creating, reading, validating, and converting BED files using bedtools (CLI) and pybedtools (Python).

## Prerequisites
```bash
# bedtools (CLI)
conda install -c bioconda bedtools

# pybedtools (Python)
pip install pybedtools
```

## Quick Start
Tell your AI agent what you want to do:
- "Create a BED file from my list of peak coordinates"
- "Validate my BED file format and fix any issues"
- "Convert my GFF file to BED format"

## Example Prompts

### Creating BED Files
> "Create a BED file from this DataFrame with columns chr, start, end, name"
> "Generate a BED4 file from my peak coordinates stored in peaks.csv"

### Validation
> "Check if my peaks.bed file has valid BED format"
> "Validate the coordinate system in my annotation.bed file"

### Conversion
> "Convert my GTF annotations to BED format"
> "Convert my VCF variants to BED intervals, adjusting for coordinate systems"

### Filtering and Sorting
> "Sort my BED file by chromosome and position"
> "Filter my BED file to keep only intervals larger than 500bp"

## What the Agent Will Do
1. Load the BED file or create one from your data
2. Validate format (check coordinates, column count, sorting)
3. Apply requested operations (filter, sort, convert)
4. Save the output with proper formatting

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

## Tips
- Always sort BED files before operations requiring sorted input
- Check chromosome naming consistency (chr1 vs 1)
- Use `pybedtools.cleanup()` after processing to remove temp files
- Validate coordinate systems when combining data from different sources
- BED files are tab-separated; spaces will cause parsing errors

## Resources
- [UCSC BED Format Specification](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)
- [bedtools Documentation](https://bedtools.readthedocs.io/)
- [pybedtools Documentation](https://daler.github.io/pybedtools/)

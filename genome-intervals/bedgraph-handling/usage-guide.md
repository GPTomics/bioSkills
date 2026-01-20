# bedGraph Handling - Usage Guide

## Overview
bedGraph is a simple text format for representing continuous genomic data like coverage or signal intensity. This skill covers creating, manipulating, and converting bedGraph files for genome browser visualization.

## Prerequisites
```bash
# bedtools
conda install -c bioconda bedtools

# UCSC tools
conda install -c bioconda ucsc-bedgraphtobigwig ucsc-bigwigtobedgraph

# pyBigWig
pip install pyBigWig

# deepTools
conda install -c bioconda deeptools
```

## Quick Start
Tell your AI agent what you want to do:
- "Generate a coverage bedGraph from my BAM file"
- "Normalize my bedGraph by library size"
- "Convert my bedGraph to bigWig format"

## Example Prompts

### Creating bedGraph Files
> "Generate a bedGraph coverage track from alignments.bam"
> "Create a normalized bedGraph from my BAM file"
> "Convert my bigWig file back to bedGraph format"

### Normalization
> "Normalize my bedGraph to reads per million"
> "Scale my coverage track by a factor of 0.5"

### Merging and Comparison
> "Merge bedGraph files from multiple samples into one file"
> "Calculate the ratio between treatment and control bedGraph files"
> "Compute the difference between two coverage tracks"

### Conversion
> "Convert my bedGraph to bigWig for genome browser"
> "Sort my bedGraph file for bigWig conversion"

### Troubleshooting
> "Fix my bedGraph file that has coordinates beyond chromosome ends"
> "Sort my bedGraph file correctly for bigWig conversion"

## What the Agent Will Do
1. Generate or load bedGraph data
2. Apply normalization or transformations if requested
3. Sort the bedGraph by chromosome and position
4. Convert to bigWig if requested
5. Verify output format is valid

## bedGraph Format
```
chr1    0       100     1.5
chr1    100     200     2.3
```
- Tab-separated
- Columns: chrom, start, end, value
- 0-based, half-open coordinates
- Must be sorted for bigWig conversion

## Key Tools

| Tool | Purpose |
|------|---------|
| bedtools genomecov | BAM to bedGraph |
| bedGraphToBigWig | bedGraph to bigWig |
| bigWigToBedGraph | bigWig to bedGraph |
| bamCoverage | BAM to normalized bigWig |

## Tips
- Always sort bedGraph before converting to bigWig: `sort -k1,1 -k2,2n`
- Use `bedClip` to fix coordinates that extend beyond chromosome ends
- Create chrom.sizes from FASTA index: `cut -f1,2 reference.fa.fai > chrom.sizes`
- Use `bedtools unionbedg` to merge multiple samples into a matrix
- bedGraph files can become very large; convert to bigWig for storage

## Resources
- [bedGraph Format](https://genome.ucsc.edu/goldenPath/help/bedgraph.html)
- [bedtools genomecov](https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html)
- [UCSC Tools](http://hgdownload.soe.ucsc.edu/admin/exe/)

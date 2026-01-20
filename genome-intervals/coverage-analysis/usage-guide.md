# Coverage Analysis - Usage Guide

## Overview
Coverage analysis calculates how well genomic regions are covered by sequencing reads. This skill covers generating bedGraph files with bedtools genomecov and calculating per-feature coverage statistics with bedtools coverage.

## Prerequisites
```bash
# bedtools
conda install -c bioconda bedtools

# pybedtools (optional)
pip install pybedtools

# samtools (for BAM operations)
conda install -c bioconda samtools
```

## Quick Start
Tell your AI agent what you want to do:
- "Generate a coverage track from my BAM file for genome browser"
- "Calculate mean coverage across my target regions"
- "Find exons with less than 20x coverage"

## Example Prompts

### Generating Coverage Tracks
> "Create a bedGraph coverage track from my alignments.bam"
> "Generate a normalized coverage track for my ChIP-seq BAM"
> "Make a coverage track scaled to reads per million"

### Per-Region Coverage
> "Calculate coverage statistics for each exon in my targets.bed"
> "What fraction of each target region is covered by reads?"
> "Summarize coverage depth across my capture regions"

### QC and Filtering
> "Find all target regions with less than 80% coverage"
> "Identify low-coverage exons in my exome sequencing data"
> "Generate a QC report of coverage across my target panel"

### Comparison
> "Compare coverage between my tumor and normal samples"
> "Generate normalized tracks for comparing ChIP-seq replicates"

## What the Agent Will Do
1. Verify BAM file is sorted and indexed
2. Generate coverage data using bedtools genomecov or coverage
3. Apply normalization if requested
4. Filter or summarize results based on thresholds
5. Output in requested format (bedGraph, summary table)

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

## Tips
- BAM files must be sorted; use `samtools sort` if needed
- Use `-bg` flag with genomecov to output bedGraph format
- Use `-scale` flag to normalize coverage (e.g., to RPM)
- For spliced alignments (RNA-seq), use `-split` flag
- Convert large bedGraph files to bigWig for efficient browser viewing
- BED input to genomecov requires a genome file (`-g` flag)

## Resources
- [bedtools genomecov](https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html)
- [bedtools coverage](https://bedtools.readthedocs.io/en/latest/content/tools/coverage.html)
- [bedGraph format](https://genome.ucsc.edu/goldenPath/help/bedgraph.html)

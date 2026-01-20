# BigWig Tracks - Usage Guide

## Overview
BigWig is an indexed binary format for continuous genomic data like coverage, ChIP-seq signal, or conservation scores. It is efficient for genome browsers and programmatic access. This skill covers creating bigWig files from bedGraph and extracting values with pyBigWig.

## Prerequisites
```bash
# UCSC tools (CLI conversion)
conda install -c bioconda ucsc-bedgraphtobigwig ucsc-bigwigtobedgraph

# pyBigWig (Python read/write)
pip install pyBigWig

# deepTools (advanced bigWig operations)
conda install -c bioconda deeptools
```

## Quick Start
Tell your AI agent what you want to do:
- "Convert my bedGraph coverage to bigWig for the genome browser"
- "Extract signal values from my bigWig file for specific regions"
- "Create a normalized bigWig track from my BAM file"

## Example Prompts

### Creating BigWig Files
> "Convert my coverage.bedGraph to bigWig format"
> "Create a normalized bigWig track from my ChIP-seq BAM"
> "Generate a CPM-normalized bigWig from alignments.bam"

### Extracting Signal
> "Get the mean signal from my bigWig file for each peak in peaks.bed"
> "Extract coverage values for the region chr1:1000000-1100000"
> "Calculate average signal across all promoters from my bigWig"

### Comparison and Analysis
> "Compare signal between two bigWig files across my regions of interest"
> "Generate a matrix of signal values for heatmap visualization"
> "Calculate fold change between treatment and control bigWig files"

### Format Conversion
> "Convert my bigWig back to bedGraph format"
> "Create a bigWig from my BAM file using deepTools"

## What the Agent Will Do
1. Prepare chromosome sizes file if needed
2. Sort bedGraph if converting from bedGraph
3. Create bigWig using appropriate tool
4. Verify output is valid and indexed
5. Extract or summarize signal if requested

## Key Concepts

### Why BigWig over bedGraph?

| Feature | bedGraph | bigWig |
|---------|----------|--------|
| File size | Large | ~10x smaller |
| Random access | Sequential only | Indexed |
| Browser support | Basic | Full |
| Query speed | Slow | Fast |

## Tips
- bedGraph must be sorted by chromosome and position before conversion
- Chromosome names in bedGraph must match exactly with chrom.sizes file
- Use `bamCoverage` from deepTools for direct BAM to normalized bigWig
- pyBigWig can query specific regions without loading the entire file
- Always provide chromosome sizes when creating bigWig files
- Use `bigWigInfo` to verify bigWig file properties

## Resources
- [pyBigWig GitHub](https://github.com/deeptools/pyBigWig)
- [UCSC Tools](http://hgdownload.soe.ucsc.edu/admin/exe/)
- [deepTools Documentation](https://deeptools.readthedocs.io/)
- [bigWig Format](https://genome.ucsc.edu/goldenPath/help/bigWig.html)

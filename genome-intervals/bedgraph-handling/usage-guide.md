# bedGraph Handling Usage Guide

This guide covers working with bedGraph files for genome browser visualization.

## bedGraph Format

```
chr1    0       100     1.5
chr1    100     200     2.3
```

- Tab-separated
- Columns: chrom, start, end, value
- 0-based, half-open coordinates
- Must be sorted for bigWig conversion

## Common Workflows

### BAM to Browser Track
1. Generate bedGraph with `bedtools genomecov`
2. Normalize by library size
3. Sort with `sort -k1,1 -k2,2n`
4. Convert to bigWig with `bedGraphToBigWig`

### Compare Samples
1. Generate normalized bedGraph for each
2. Merge with `bedtools unionbedg`
3. Calculate ratios or differences

## Requirements

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

## Key Tools

| Tool | Purpose |
|------|---------|
| bedtools genomecov | BAM to bedGraph |
| bedGraphToBigWig | bedGraph to bigWig |
| bigWigToBedGraph | bigWig to bedGraph |
| bamCoverage | BAM to normalized bigWig |

## Troubleshooting

### "bedGraph not sorted"
Run: `sort -k1,1 -k2,2n input.bedgraph > sorted.bedgraph`

### "bedGraph extends beyond chromosome"
Run: `bedClip input.bedgraph chrom.sizes clipped.bedgraph`

### Missing chromosomes
Generate chrom.sizes: `cut -f1,2 reference.fa.fai > chrom.sizes`

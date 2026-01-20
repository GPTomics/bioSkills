# Amplicon Processing - Usage Guide

## Overview

DADA2 infers exact amplicon sequence variants (ASVs) from amplicon sequencing data, providing single-nucleotide resolution without clustering.

## Prerequisites

```bash
# R packages
install.packages('BiocManager')
BiocManager::install('dada2')
```

## Quick Start

Tell your AI agent what you want to do:
- "Process my 16S amplicon data with DADA2"
- "Denoise paired-end reads and generate an ASV table"

## Example Prompts

### Basic Processing
> "Run DADA2 on my paired-end 16S data in /data/fastq/"

> "Generate quality plots for my amplicon sequences"

### Parameter Tuning
> "What truncLen should I use based on these quality profiles?"

> "Process my reads with maxEE=2 and truncLen of 240,160"

### Troubleshooting
> "My merge rate is low, help me diagnose the issue"

> "Remove chimeras from my ASV table"

## What the Agent Will Do

1. Inspect quality profiles to determine trim points
2. Filter and trim reads to remove low-quality bases
3. Learn error rates from the data
4. Denoise reads to infer true sequences
5. Merge paired reads (if applicable)
6. Remove chimeric sequences
7. Generate final ASV table and sequences

## Tips

- ASVs provide single-nucleotide resolution vs 97% OTU clustering
- Examine quality profiles before choosing truncLen
- Ensure sufficient overlap for merging (20bp minimum)
- Normal chimera rate is 10-25% of ASVs but few reads
- Save output as RDS for downstream phyloseq analysis

## ASVs vs OTUs

| Feature | ASVs (DADA2) | OTUs (97%) |
|---------|--------------|------------|
| Resolution | Single nucleotide | 3% similarity |
| Reproducibility | Exact sequences | Cluster-dependent |
| Comparability | Can merge studies | Requires re-clustering |
| Sensitivity | Higher | Lower |

## Key Parameters

### filterAndTrim

| Parameter | Description | Typical Value |
|-----------|-------------|---------------|
| truncLen | Truncate reads at position | Based on quality plots |
| maxEE | Maximum expected errors | 2 |
| maxN | Maximum Ns allowed | 0 |
| truncQ | Truncate at quality score | 2 |

## Common Issues

### Low merge rate
- Check amplicon length vs read length
- Increase minimum overlap
- Reads may not overlap (too long amplicon)

### High chimera rate
- Normal: 10-25% of ASVs, but few reads
- Check PCR cycles (lower is better)

### Few reads passing filter
- Loosen maxEE (try 2, 5)
- Adjust truncLen
- Check raw data quality

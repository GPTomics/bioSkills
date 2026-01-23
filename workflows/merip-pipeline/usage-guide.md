# MeRIP-seq Pipeline - Usage Guide

## Overview

Complete workflow from MeRIP-seq FASTQ to m6A peaks, differential methylation, and visualization.

## Prerequisites

```bash
conda install -c bioconda star samtools
```

```r
BiocManager::install(c('exomePeak2', 'Guitar', 'TxDb.Hsapiens.UCSC.hg38.knownGene'))
```

## Quick Start

- "Analyze my MeRIP-seq data for m6A peaks"
- "Run the m6A pipeline from FASTQ to differential methylation"
- "Process my epitranscriptomics data end-to-end"

## Example Prompts

### Full Pipeline

> "Run the complete MeRIP-seq pipeline"

> "Find m6A peaks and compare between conditions"

### Specific Steps

> "Just call peaks from my aligned BAMs"

> "Create metagene plot of m6A distribution"

## What the Agent Will Do

1. Align IP and Input samples (STAR)
2. QC alignment and IP enrichment
3. Call m6A peaks (exomePeak2)
4. Run differential methylation if multiple conditions
5. Generate metagene and motif plots
6. Export peaks and statistics

## Tips

- **Paired samples** - Match IP with corresponding Input
- **Replicates** - At least 2 biological replicates per condition
- **DRACH motif** - m6A consensus; peaks should be enriched
- **Stop codon** - Classic m6A enrichment pattern
- **FDR < 0.05** - Standard significance threshold

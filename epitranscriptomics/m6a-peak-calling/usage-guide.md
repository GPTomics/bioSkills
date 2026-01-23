# m6A Peak Calling - Usage Guide

## Overview

Identify m6A modification sites by comparing MeRIP-seq IP enrichment to input background.

## Prerequisites

```r
BiocManager::install(c('exomePeak2', 'MeTPeak'))
```

```bash
conda install -c bioconda macs3
```

## Quick Start

- "Call m6A peaks from my MeRIP-seq data"
- "Find m6A sites using exomePeak2"
- "Identify significant IP enrichment regions"

## Example Prompts

### Basic Peak Calling

> "Run exomePeak2 on my IP and input BAM files"

> "Call m6A peaks with MACS3"

### With Replicates

> "Call peaks with 3 IP and 3 input replicates"

> "Identify m6A sites requiring concordance across replicates"

## What the Agent Will Do

1. Load IP and input BAM files
2. Perform peak calling (exomePeak2 or MACS3)
3. Apply significance thresholds
4. Export peaks in BED format
5. Generate peak statistics

## Tips

- **exomePeak2** - Recommended for transcript-aware analysis
- **MACS3** - Faster, good for initial exploration
- **Replicates** - Use at least 2 per condition
- **FDR < 0.05** - Standard significance threshold
- **DRACH motif** - m6A occurs at DRACH consensus (D=A/G/U, R=A/G, H=A/C/U)

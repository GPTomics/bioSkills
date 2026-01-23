# MeRIP-seq Preprocessing - Usage Guide

## Overview

Align MeRIP-seq IP and input samples and perform QC for m6A peak calling.

## Prerequisites

```bash
conda install -c bioconda star samtools deeptools
```

## Quick Start

- "Align my MeRIP-seq samples"
- "Check IP enrichment quality"
- "Assess replicate correlation for MeRIP"

## Example Prompts

### Alignment

> "Align my MeRIP-seq IP and input FASTQ files to the human genome"

> "Process paired-end MeRIP-seq data with STAR"

### Quality Control

> "Check correlation between my MeRIP replicates"

> "Assess IP enrichment in my m6A experiment"

## What the Agent Will Do

1. Build STAR index if needed
2. Align IP and input samples
3. Sort and index BAM files
4. Check replicate correlation
5. Assess IP/input signal patterns

## Tips

- **Paired samples** - IP and input from same biological replicate
- **Replicates** - At least 2 biological replicates per condition
- **Sequencing depth** - 30-50M reads per sample typical
- **Strand-specific** - Most MeRIP protocols are unstranded

# CLIP-seq Alignment - Usage Guide

## Overview

Align preprocessed CLIP-seq reads to the genome optimized for crosslink site identification.

## Prerequisites

```bash
conda install -c bioconda star bowtie2 samtools
```

## Quick Start

- "Align my CLIP reads with STAR"
- "Use unique mapping only"
- "Deduplicate after alignment"

## Example Prompts

> "Align CLIP reads allowing 1 mismatch"

> "Map with bowtie2 in very-sensitive mode"

## What the Agent Will Do

1. Align reads with STAR or bowtie2
2. Filter for unique mappers
3. Sort and index BAM
4. Deduplicate with UMI tools

## Tips

- **Unique mapping** preferred for precise binding sites
- **EndToEnd** alignment prevents soft-clipping

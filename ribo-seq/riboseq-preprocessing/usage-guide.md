# Ribo-seq Preprocessing - Usage Guide

## Overview

Preprocess ribosome profiling data including adapter trimming, size selection for ribosome-protected fragments, rRNA depletion, and alignment.

## Prerequisites

```bash
conda install -c bioconda cutadapt sortmerna bowtie2 star samtools
```

## Quick Start

Tell your AI agent:
- "Preprocess my Ribo-seq FASTQ files"
- "Remove rRNA contamination from my reads"
- "Filter reads to ribosome footprint length"
- "Align Ribo-seq reads to the transcriptome"

## Example Prompts

### Trimming and Filtering

> "Trim adapters and select 28-32 nt reads from my Ribo-seq data"

> "What percentage of reads are in the footprint size range?"

> "Filter my trimmed reads by length"

### rRNA Removal

> "Remove rRNA contamination with SortMeRNA"

> "What fraction of reads map to rRNA?"

> "Build an rRNA index for depletion"

### Alignment

> "Align Ribo-seq reads to the human transcriptome"

> "Generate a BAM file with uniquely mapped reads only"

> "Check alignment statistics"

## What the Agent Will Do

1. Trim 3' adapters from raw reads
2. Filter to ribosome footprint size (28-32 nt)
3. Remove rRNA-mapping reads
4. Align remaining reads to transcriptome
5. Generate sorted, indexed BAM

## Tips

- **28-32 nt reads** are typical ribosome footprints; adjust per organism
- **rRNA removal** is critical - Ribo-seq has high rRNA contamination
- **Unique mapping** recommended for precise positioning
- **Read length peak** should be ~30 nt for good library
- **Multimapping** - disable or handle carefully for paralogs

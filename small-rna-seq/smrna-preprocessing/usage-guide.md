# Small RNA Preprocessing - Usage Guide

## Overview

Preprocess small RNA-seq data with adapter trimming and size selection optimized for miRNA, piRNA, and snoRNA analysis.

## Prerequisites

```bash
pip install cutadapt
# or
conda install -c bioconda cutadapt fastp
```

## Quick Start

Tell your AI agent:
- "Trim adapters from my small RNA-seq data"
- "Filter reads to miRNA length range (18-26 nt)"
- "Collapse identical sequences for miRNA quantification"
- "Check read length distribution in my small RNA library"

## Example Prompts

### Adapter Trimming

> "Trim Illumina TruSeq small RNA adapters from my FASTQ files"

> "Remove adapters and filter to 18-30 nt length"

> "Process my small RNA samples with quality trimming and adapter removal"

### Size Selection

> "Filter my trimmed reads to keep only miRNA-sized fragments (18-26 nt)"

> "Separate piRNA-sized reads (26-32 nt) from miRNAs"

> "What percentage of my reads are in the miRNA size range?"

### Quality Control

> "Plot the read length distribution after trimming"

> "Show me adapter trimming statistics"

> "Check if my small RNA library is good quality"

## What the Agent Will Do

1. Identify adapter sequence (from kit or auto-detect)
2. Trim 3' adapters with cutadapt or fastp
3. Apply size selection filters (18-26 nt for miRNA)
4. Optionally collapse identical sequences
5. Generate QC report with length distribution

## Tips

- **Specify adapter sequence** - auto-detection may miss partial adapters
- **miRNA peak at 21-23 nt** - if peak is elsewhere, check library prep
- **High adapter content is good** - indicates successful size selection
- **Collapse reads** for faster downstream analysis
- **Keep quality > 20** to reduce noise in quantification

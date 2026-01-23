# MiXCR Analysis - Usage Guide

## Overview

Align and assemble immune receptor sequences to identify TCR and BCR clonotypes from bulk or single-cell sequencing data.

## Prerequisites

```bash
conda install -c bioconda mixcr
# or download from https://github.com/milaboratory/mixcr/releases
```

## Quick Start

Tell your AI agent:
- "Process my TCR-seq FASTQ files with MiXCR"
- "Assemble clonotypes from immune repertoire data"
- "Export clonotypes with CDR3 sequences"
- "Run MiXCR on my 10x VDJ data"

## Example Prompts

### Complete Workflow

> "Run MiXCR amplicon analysis on my paired-end TCR FASTQ files"

> "Process BCR sequencing data and export clonotypes"

> "Analyze 10x Genomics VDJ data with MiXCR"

### Custom Analysis

> "Align my repertoire reads to human V(D)J genes"

> "Assemble clonotypes allowing partial alignments"

> "Export clones with V, D, J gene usage information"

### Quality Control

> "Show me MiXCR alignment statistics"

> "What percentage of reads have CDR3 identified?"

> "Check clonotype assembly quality"

## What the Agent Will Do

1. Select appropriate MiXCR preset based on data type
2. Run alignment to V(D)J reference
3. Assemble partial alignments if needed
4. Assemble clonotypes
5. Export results with CDR3 sequences and gene usage

## Tips

- **Use presets** for common protocols (amplicon, 10x-vdj-tcr)
- **Check alignment rate** - >80% is good, <50% suggests problems
- **CDR3 is key** - most analyses focus on CDR3 diversity
- **Paired-end recommended** for full-length CDR3 assembly
- **Memory** - large datasets need 16GB+ RAM

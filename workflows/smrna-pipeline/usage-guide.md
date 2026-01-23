# Small RNA-seq Pipeline - Usage Guide

## Overview

Complete workflow from small RNA-seq FASTQ to differential miRNA expression and target prediction.

## Prerequisites

```bash
# CLI tools
conda install -c bioconda cutadapt mirdeep2 miranda

# R packages
BiocManager::install('DESeq2')
```

## Quick Start

- "Analyze my small RNA-seq data for differential miRNAs"
- "Run the miRNA pipeline from FASTQ to targets"
- "Process my miRNA sequencing end-to-end"

## Example Prompts

### Full Pipeline

> "Run the complete small RNA-seq pipeline"

> "Find differentially expressed miRNAs between my conditions"

### Specific Steps

> "Just run miRDeep2 quantification"

> "Predict targets for my differentially expressed miRNAs"

## What the Agent Will Do

1. Trim adapters and filter by size (cutadapt)
2. Align to genome (mapper.pl)
3. Quantify known and novel miRNAs (miRDeep2)
4. Run differential expression (DESeq2)
5. Predict mRNA targets (miRanda)
6. Generate QC and summary reports

## Tips

- **Adapter** - Check your kit's 3' adapter sequence
- **Size filter** - 18-30nt captures miRNAs (21-23nt peak)
- **Novel miRNAs** - Require structural validation
- **Target prediction** - Filter by binding energy and conservation
- **Multiple testing** - Use adjusted p-values (padj < 0.05)

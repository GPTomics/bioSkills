# miRDeep2 Analysis - Usage Guide

## Overview

Discover novel miRNAs and quantify known miRNAs using miRDeep2's de novo prediction algorithm based on secondary structure and read patterns.

## Prerequisites

```bash
conda install -c bioconda mirdeep2 bowtie

# Download miRBase references
wget https://www.mirbase.org/download/mature.fa
wget https://www.mirbase.org/download/hairpin.fa
```

## Quick Start

Tell your AI agent:
- "Discover novel miRNAs in my small RNA-seq data"
- "Run miRDeep2 on my collapsed reads"
- "Quantify known human miRNAs with miRDeep2"
- "Set up miRDeep2 with the human genome"

## Example Prompts

### Novel miRNA Discovery

> "Find novel miRNAs in my collapsed FASTA using miRDeep2"

> "Run miRDeep2 discovery with human miRBase references"

> "Predict novel miRNAs and score them by confidence"

### Quantification

> "Quantify known human miRNAs in my samples"

> "Run miRDeep2 quantifier on multiple samples"

> "Get expression counts for mature miRNAs"

### Results Analysis

> "Parse miRDeep2 results and filter high-confidence predictions"

> "How many novel miRNAs have score > 10?"

> "Compare my novel miRNAs to miRBase"

## What the Agent Will Do

1. Build bowtie index for genome (if needed)
2. Run mapper.pl to align collapsed reads
3. Execute miRDeep2.pl for discovery and quantification
4. Parse output files for novel miRNA predictions
5. Filter by confidence score

## Tips

- **Collapse reads first** - miRDeep2 requires pre-collapsed FASTA
- **miRDeep2 score > 10** is high confidence for novel miRNAs
- **Check randfold p-value** - low p indicates stable hairpin structure
- **Use species-specific miRBase** - improves known miRNA mapping
- **Run time** - can take hours for large datasets; use multiple cores

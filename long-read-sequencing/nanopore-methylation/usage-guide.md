# Nanopore Methylation Calling Usage Guide

## Overview

This guide covers detecting DNA methylation directly from Oxford Nanopore sequencing data.

## Prerequisites

```bash
# modkit (ONT tool)
# Download from https://github.com/nanoporetech/modkit

# Dorado basecaller with methylation models
# Download from https://github.com/nanoporetech/dorado

# Supporting tools
conda install -c bioconda minimap2 samtools
```

## Quick Start

Tell your AI agent what you want to do:
- "Call 5mC methylation from my nanopore BAM file"
- "Generate a methylation bedMethyl file from my ONT data"
- "Compare methylation between tumor and normal samples"
- "Extract methylation in promoter regions"

## Example Prompts

### Basic Methylation Calling

> "Run modkit pileup on my aligned BAM to get CpG methylation calls"

> "Generate a bedMethyl file combining both strands for CpG sites"

### Region Analysis

> "Extract methylation levels in CpG islands from my nanopore data"

> "Calculate average methylation in gene promoters"

### Differential Analysis

> "Compare methylation between my two samples and find differentially methylated regions"

## What the Agent Will Do

1. Verify BAM contains methylation tags (MM/ML)
2. Run modkit pileup with appropriate settings
3. Filter by coverage and quality thresholds
4. Generate bedMethyl output
5. Summarize methylation statistics

## Tips

- Ensure basecalling used methylation-aware model
- Minimum 10x coverage recommended for reliable calls
- Combine strands for CpG (--combine-strands)
- Use --cpg flag to limit to CpG context
- Check modkit summary before analysis
- For 6mA detection, use appropriate model and flags

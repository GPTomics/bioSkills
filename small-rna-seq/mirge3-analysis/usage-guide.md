# miRge3 Analysis - Usage Guide

## Overview

Fast miRNA quantification with miRge3, including isomiR detection, A-to-I RNA editing analysis, and multi-species support.

## Prerequisites

```bash
pip install mirge3

# Download organism library
miRge3.0 --download-library human mirbase
```

## Quick Start

Tell your AI agent:
- "Quantify miRNAs with miRge3"
- "Detect isomiRs in my small RNA-seq"
- "Run miRge3 on multiple samples"
- "Analyze A-to-I editing in my miRNAs"

## Example Prompts

### Basic Quantification

> "Run miRge3 on my trimmed FASTQ files"

> "Quantify human miRNAs using miRBase"

> "Process all my samples and create a count matrix"

### IsomiR Analysis

> "Detect isomiR variants for each miRNA"

> "Find the dominant isomiR for each miRNA family"

> "Analyze 3' non-templated additions"

### RNA Editing

> "Detect A-to-I editing sites in my miRNAs"

> "What percentage of reads show editing at known sites?"

> "Compare editing levels between samples"

## What the Agent Will Do

1. Check for miRge3 library installation
2. Configure analysis parameters (organism, adapter)
3. Run miRge3 annotate on input FASTQs
4. Parse count matrices and QC report
5. Optionally analyze isomiRs or editing

## Tips

- **miRge3 is faster** than miRDeep2 for known miRNA quantification
- **Use --isomir** flag to get isomiR-level resolution
- **Libraries must match** your reference genome version
- **Adapter trimming** can be done internally or externally
- **RPM normalization** provided automatically for comparison

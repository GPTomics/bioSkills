# ORF Detection - Usage Guide

## Overview

Detect actively translated open reading frames from Ribo-seq data, including canonical ORFs, upstream ORFs (uORFs), and novel translated regions.

## Prerequisites

```bash
pip install RiboCode
# or
conda install -c bioconda ribocode riborf

# For ORF quantification (R)
# BiocManager::install('ORFik')
```

## Quick Start

Tell your AI agent:
- "Find translated ORFs in my Ribo-seq data"
- "Detect uORFs in my genes"
- "Run RiboCode to identify novel translated regions"
- "Quantify ORF-level translation with ORFquant"

## Example Prompts

### ORF Detection

> "Run RiboCode on my Ribo-seq BAM file"

> "Find all translated ORFs with p-value < 0.05"

> "How many novel ORFs are detected?"

### uORF Analysis

> "Find upstream ORFs in my genes of interest"

> "Which uORFs are actively translated?"

> "Compare uORF translation between conditions"

### Classification

> "Classify detected ORFs by type (annotated, uORF, novel)"

> "Find translated regions in long non-coding RNAs"

### Quantification

> "Quantify ORF-level translation with ORFik"

> "Compare ORF expression across conditions"

## What the Agent Will Do

1. Prepare annotation files for RiboCode
2. Run ORF detection with P-site offset correction
3. Filter by statistical significance
4. Classify ORFs by genomic context
5. Report translated ORFs with coordinates

## Tips

- **P-site offset** must be correct for accurate detection
- **Periodicity** in detected ORFs validates true translation
- **uORFs** can regulate downstream CDS translation
- **Novel ORFs** may encode micropeptides
- **Binomial test** assesses 3-nt periodicity significance
- **ORFquant/ORFik** provides quantitative ORF expression for differential analysis

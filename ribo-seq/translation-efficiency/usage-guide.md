# Translation Efficiency - Usage Guide

## Overview

Calculate translation efficiency (TE) as the ratio of ribosome occupancy to mRNA abundance, enabling detection of translational regulation independent of transcription.

## Prerequisites

```r
# R
BiocManager::install(c('riborex', 'DESeq2'))
```

```bash
# Python
pip install plastid pandas numpy
```

## Quick Start

Tell your AI agent:
- "Calculate translation efficiency from my Ribo-seq and RNA-seq"
- "Find genes with differential TE between conditions"
- "Run riborex on my count matrices"
- "Identify translationally regulated genes"

## Example Prompts

### Calculate TE

> "Calculate TE for each gene using Ribo-seq and RNA-seq counts"

> "Normalize counts and compute log2 translation efficiency"

> "Create a TE matrix for all samples"

### Differential TE

> "Find genes with significantly different TE between treatment and control"

> "Run riborex to detect differential translation"

> "Use DESeq2 interaction model for TE analysis"

### Interpretation

> "Which genes are translationally upregulated?"

> "Find genes with increased mRNA but decreased TE"

> "Identify genes regulated only at translation level"

## What the Agent Will Do

1. Load Ribo-seq and RNA-seq count matrices
2. Normalize counts (TPM or RPKM)
3. Calculate TE per gene per sample
4. Run statistical test for differential TE
5. Filter by significance (padj < 0.05)

## Tips

- **Match samples** - Ribo-seq and RNA-seq must be from same biological samples
- **CDS counts for Ribo-seq** - only count in coding sequence
- **Full transcript for RNA-seq** - include UTRs
- **Pseudocounts** avoid division by zero
- **Log2 TE** is more interpretable for fold changes

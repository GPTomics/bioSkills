# CLIP Peak Calling - Usage Guide

## Overview

Call protein-RNA binding site peaks from aligned CLIP-seq data.

## Prerequisites

```bash
conda install -c bioconda clipper piranha pureclip
```

## Quick Start

Tell your AI agent what you want to do:
- "Call peaks with CLIPper"
- "Run PureCLIP for HMM-based crosslink site detection"
- "Find significant binding clusters"
- "Call peaks with input control normalization"

## Example Prompts

### CLIPper

> "Run CLIPper on my deduped BAM"

> "Call peaks with FDR < 0.05"

### PureCLIP

> "Run PureCLIP for single-nucleotide crosslink detection"

> "Call peaks with PureCLIP using my input control"

> "Get high-confidence binding regions from eCLIP data"

## What the Agent Will Do

1. Run peak caller (CLIPper, PureCLIP, or Piranha) on aligned BAM
2. Apply significance filtering (FDR, p-value, or score threshold)
3. Output BED files with peaks or crosslink sites
4. Calculate quality metrics (FRiP, peak count)

## Tips

- **CLIPper** is ENCODE standard for eCLIP
- **PureCLIP** uses HMM for single-nucleotide crosslink site resolution
- **Use deduplicated** BAM for accurate calling
- **Input control** improves specificity (use with PureCLIP -ibam)
- **Score filtering**: PureCLIP score >= 3 for medium confidence, >= 5 for high

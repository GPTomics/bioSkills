# CLIP Peak Calling - Usage Guide

## Overview

Call protein-RNA binding site peaks from aligned CLIP-seq data.

## Prerequisites

```bash
conda install -c bioconda clipper piranha
```

## Quick Start

- "Call peaks with CLIPper"
- "Find significant binding clusters"

## Example Prompts

> "Run CLIPper on my deduped BAM"

> "Call peaks with p-value < 0.01"

## What the Agent Will Do

1. Run peak caller on aligned BAM
2. Filter by significance
3. Output BED file with peaks

## Tips

- **CLIPper** is ENCODE standard for eCLIP
- **Use deduplicated** BAM for accurate calling

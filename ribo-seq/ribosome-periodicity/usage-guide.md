# Ribosome Periodicity - Usage Guide

## Overview

Validate Ribo-seq data quality by checking 3-nucleotide periodicity (the hallmark of active translation) and determining P-site offsets for accurate ribosome positioning.

## Prerequisites

```bash
pip install plastid numpy matplotlib scipy
```

## Quick Start

Tell your AI agent:
- "Check 3-nucleotide periodicity in my Ribo-seq data"
- "Determine P-site offset for my read lengths"
- "Create a metagene plot around start codons"
- "Validate my Ribo-seq library quality"

## Example Prompts

### Periodicity Validation

> "Check if my Ribo-seq data shows 3-nt periodicity"

> "Is my ribosome profiling library good quality?"

> "Plot periodicity around start and stop codons"

### P-site Offset

> "Calculate optimal P-site offset for 28-30 nt reads"

> "What offset should I use for my read length distribution?"

> "Generate offset lookup table by read length"

### Quality Assessment

> "Create a metagene plot for my Ribo-seq"

> "Compare periodicity between samples"

> "What fraction of reads show clear periodicity?"

## What the Agent Will Do

1. Load aligned Ribo-seq BAM file
2. Extract reads around annotated start codons
3. Calculate read density by position and frame
4. Quantify 3-nt periodicity (FFT or frame enrichment)
5. Determine optimal P-site offset per read length

## Tips

- **Strong periodicity** in frame 0 indicates good library quality
- **P-site offset** varies by read length - calculate separately
- **12-13 nt offset** is typical for 28-30 nt footprints
- **Weak periodicity** may indicate poor digestion or library issues
- **Use annotated CDS** for metagene analysis

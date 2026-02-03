# qPCR Primer and Probe Design - Usage Guide

## Overview

This skill covers designing primers and probes for quantitative PCR (real-time PCR) using primer3-py. Supports TaqMan hydrolysis probes and SYBR Green primer-only assays.

## Prerequisites

```bash
pip install primer3-py biopython pandas
```

## Quick Start

Ask your AI agent:

- "Design qPCR primers with a TaqMan probe for my gene"
- "Design SYBR Green primers for this sequence"
- "Find qPCR primers that span exon 2-3 junction"

## Example Prompts

### TaqMan Assay
> "Design primers and a TaqMan probe for quantifying GAPDH"

> "Create a qPCR assay with probe Tm around 70C"

### SYBR Green Assay
> "Design SYBR Green primers for a 100bp amplicon"

> "Find qPCR primers with minimal dimer potential"

### cDNA-Specific
> "Design primers that span the exon junction so they won't amplify genomic DNA"

> "Create exon-spanning primers for my mRNA target"

### Multiplex
> "Design qPCR primers for both my target gene and housekeeping gene with matching Tms"

## What the Agent Will Do

1. Load your template sequence
2. Configure primer3 for qPCR-specific parameters
3. Design primers with optimal Tm (~60C)
4. Design probe with elevated Tm (~70C) if requested
5. Return primer/probe sets ranked by quality

## Tips

- **Short amplicons** - 70-150bp is ideal for qPCR efficiency
- **Probe Tm** - Should be 8-10C higher than primer Tm
- **Avoid G at probe 5'** - G adjacent to FAM quenches fluorescence
- **Exon spanning** - Prevents genomic DNA amplification in RNA quantification
- **Tm matching** - Keep primer Tms within 1-2C for consistent annealing
- **SYBR specificity** - Use strict complementarity settings to minimize primer dimers

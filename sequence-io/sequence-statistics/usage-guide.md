# Sequence Statistics - Usage Guide

## Overview

This skill enables AI agents to help you calculate comprehensive statistics for sequence datasets using Biopython.

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:

- "Calculate N50 for my assembly"
- "Generate a summary report for these sequences"
- "Compare statistics across multiple FASTA files"
- "Show me the GC content distribution"

## Example Prompts

### Assembly QC
> "Calculate N50, L50, and total bases for my genome assembly"

### Comparison
> "Compare assembly statistics for all FASTA files in the assemblies folder"

### Distribution
> "Show me the length distribution histogram for my sequences"

### Export
> "Generate a CSV with statistics for all my sequence files"

## Key Metrics

| Metric | What It Tells You |
|--------|-------------------|
| N50 | Assembly contiguity (higher = better) |
| L50 | How many sequences contain half the data |
| Total bp | Dataset size |
| GC% | Nucleotide composition |

## What the Agent Will Do
1. Load sequences from input file
2. Calculate requested statistics (length, GC, etc.)
3. Aggregate statistics across all sequences
4. Report summary results

## Tips

- N50 is the most common assembly quality metric
- Compare N50 across different assembly parameters
- GC content varies by organism (humans ~41%, E. coli ~51%)
- Use CSV export for downstream analysis in R/Python

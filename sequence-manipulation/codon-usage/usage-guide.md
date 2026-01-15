# Codon Usage - Usage Guide

## Overview

This skill enables AI agents to help you analyze codon usage patterns in coding sequences using Biopython. It covers codon counting, CAI calculation, RSCU analysis, and codon optimization for heterologous expression.

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:

- "Calculate codon frequencies for this coding sequence"
- "What is the CAI of this gene for E. coli expression?"
- "Find rare codons in my sequence"
- "Optimize this gene for yeast expression"

## Example Prompts

### Basic Codon Analysis
> "Count the codons in this coding sequence and show frequencies"

### CAI Calculation
> "Calculate the Codon Adaptation Index for this gene using E. coli codon usage"

### RSCU Analysis
> "Calculate RSCU values for my coding sequence"

### Rare Codon Detection
> "Find codons that are used less than 10% of the time"

### Codon Optimization
> "Replace rare codons with preferred synonymous codons for E. coli"

### Compare Sequences
> "Compare codon usage between these two genes"

### Wobble Position Bias
> "What is the GC content at each codon position?"

## What the Agent Will Do

1. Import Bio.SeqUtils.CodonUsage and related modules
2. Parse the coding sequence
3. Calculate requested metrics (CAI, RSCU, frequencies)
4. Return formatted analysis results
5. Suggest optimizations if requested

## Key Concepts

### CAI (Codon Adaptation Index)
- Measures how well-adapted a gene's codon usage is to a host organism
- Range: 0-1 (higher = better adapted)
- Requires a reference set of highly expressed genes

### RSCU (Relative Synonymous Codon Usage)
- Measures bias toward specific synonymous codons
- 1.0 = no bias (all synonymous codons used equally)
- \>1 = overused, <1 = underused

### Nc (Effective Number of Codons)
- Measures overall codon bias
- Range: 20-61
- Lower values = more biased usage

## Tips

- Ensure sequences are in frame (length divisible by 3)
- Stop codons are usually excluded from analysis
- CAI requires training on reference genes from target organism
- Consider organism-specific codon tables for non-standard genetic codes
- GC3 (GC at wobble position) correlates with overall genome GC content

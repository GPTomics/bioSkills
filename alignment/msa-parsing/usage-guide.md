# MSA Parsing - Usage Guide

## Overview

This skill focuses on parsing and analyzing multiple sequence alignments (MSAs). It covers extracting information, analyzing gaps and conservation, filtering sequences, and preparing alignments for downstream analysis like phylogenetics or structure prediction.

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:
- "Find all the fully conserved positions in this alignment"
- "Remove sequences with more than 10% gaps"
- "Generate a consensus sequence from this alignment"

## Example Prompts

### Analyzing Content
> "Show me the composition of each column in the alignment"

> "Find positions that are conserved in at least 80% of sequences"

> "Count the gaps in each sequence"

### Filtering and Cleaning
> "Remove columns with more than 50% gaps"

> "Filter out sequences with too many gaps"

> "Remove duplicate sequences from the alignment"

### Extracting Information
> "Get the sequence for species_A from this alignment"

> "Extract columns 100-200 from the alignment"

> "List all sequence IDs in the alignment"

### Consensus and Conservation
> "Generate a consensus sequence with 70% threshold"

> "Find the most conserved regions in this alignment"

> "What is the consensus at each position?"

## What the Agent Will Do

1. Load the alignment file
2. Parse alignment structure (sequences, columns, gaps)
3. Perform requested analysis or filtering
4. Return results (statistics, filtered alignment, consensus)
5. Optionally save modified alignment

## Key Concepts

| Term | Description |
|------|-------------|
| Column | Vertical slice (same position across all sequences) |
| Conservation | Fraction of sequences with same residue at a position |
| Consensus | Most common character at each position |
| Gap | Missing data represented by '-' |

## Tips

- Always check gap content before phylogenetic analysis
- Gappy columns can be artifacts of alignment algorithms
- Conservation thresholds depend on alignment diversity
- Stockholm format preserves annotations other formats lose
- Position numbering is 0-based (first column is 0)

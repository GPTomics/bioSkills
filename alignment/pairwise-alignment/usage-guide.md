# Pairwise Alignment - Usage Guide

## Overview

This skill performs pairwise sequence alignment to compare two DNA, RNA, or protein sequences. It uses Biopython's `PairwiseAligner` class which implements dynamic programming algorithms for finding optimal alignments.

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:

- "Align these two DNA sequences and show me the best alignment"
- "Compare this protein sequence against a reference using BLOSUM62"
- "Find the best matching region between these two sequences"

## Example Prompts

### Global Alignment
> "Perform a global alignment between ACCGGTAACGTAG and ACCGTTAACGAAG"

> "Align the first two sequences in my FASTA file"

### Local Alignment
> "Find the best local alignment between these sequences to identify conserved regions"

> "Use Smith-Waterman to find matching regions in these proteins"

### Protein Alignment
> "Align these two protein sequences using BLOSUM62 scoring"

> "Compare my query protein against the reference with appropriate gap penalties"

### Scoring and Analysis
> "Calculate the alignment score between these sequences"

> "Show me all optimal alignments and their scores"

## What the Agent Will Do

1. Create a `PairwiseAligner` with appropriate settings
2. Configure scoring (match/mismatch for DNA, substitution matrix for proteins)
3. Set gap penalties based on sequence type
4. Generate optimal alignment(s)
5. Display aligned sequences with match indicators
6. Report alignment score and statistics

## Supported Alignment Modes

| Mode | Best For |
|------|----------|
| Global | Full-length sequence comparison |
| Local | Finding conserved regions |
| Semiglobal | Overlapping sequences, primer alignment |

## Tips

- For **DNA/RNA**: Use simple match/mismatch scoring (match=2, mismatch=-1)
- For **proteins**: Always use a substitution matrix like BLOSUM62
- For **finding regions**: Use local mode instead of global
- For **many alignments**: Set `max_alignments` to limit memory usage
- Alignment score depends heavily on gap penalties - adjust if results seem wrong

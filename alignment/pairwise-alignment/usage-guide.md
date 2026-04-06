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
| Global | Full-length comparison of similar-length homologs |
| Local | Finding conserved domains, divergent termini, database-style searches |
| Semiglobal | Fragment-to-reference, overlap detection, primer alignment |

## Tips

- For **DNA/RNA**: Use simple match/mismatch scoring (match=2, mismatch=-1); or load NUC.4.4 for IUPAC ambiguity support
- For **proteins**: Always use a substitution matrix. BLOSUM62 for general use, BLOSUM80 for close homologs, BLOSUM45 for distant
- For **finding regions**: Use local mode instead of global; global alignment of very different-length sequences forces meaningless terminal gaps
- For **many alignments**: Set `max_alignments` to limit memory usage
- Gap penalties heavily impact results: always use affine (gap open >> gap extend) to model real indel biology; typical protein values: open=-11, extend=-1
- Sequences below ~25% protein identity are in the **twilight zone** where alignment reliability drops sharply; consider profile methods (HHpred) or structural alignment
- Percent identity has multiple definitions (different denominators); always report which method was used
- Alignment algorithms always produce output even for unrelated sequences. Check statistical significance (E-value/bit score) to confirm homology

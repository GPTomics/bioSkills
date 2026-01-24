# Alignment Statistics - Usage Guide

## Overview

This skill calculates statistical metrics for sequence alignments including identity, conservation, entropy, and substitution patterns. These metrics are essential for assessing alignment quality, identifying conserved regions, and understanding evolutionary relationships.

## Prerequisites

```bash
pip install biopython numpy
```

## Quick Start

Tell your AI agent what you want to do:
- "Calculate pairwise identity between all sequences in this alignment"
- "Show me the conservation score for each column"
- "What is the average sequence identity in this alignment?"

## Example Prompts

### Identity Calculations
> "Create a pairwise identity matrix for this alignment"

> "What is the percent identity between sequence A and sequence B?"

> "Find the most similar pair of sequences in the alignment"

### Conservation Analysis
> "Calculate the conservation score at each position"

> "Which columns are most conserved?"

> "Plot a conservation profile across the alignment"

### Information Content
> "Calculate Shannon entropy for each column"

> "What is the information content at each position?"

> "Find the most variable positions in the alignment"

### Gap Analysis
> "What fraction of the alignment is gaps?"

> "Which sequences have the most gaps?"

> "How many gap-free columns are there?"

### Substitution Patterns
> "Count the substitutions between all pairs of sequences"

> "What are the most common substitution types?"

> "Build a substitution matrix from this alignment"

## What the Agent Will Do

1. Load the alignment file
2. Calculate requested metrics (identity, conservation, entropy, etc.)
3. Summarize results (averages, distributions, extremes)
4. Identify notable patterns (highly conserved/variable regions)
5. Output tables, matrices, or profiles as appropriate

## Key Metrics Explained

| Metric | What It Measures | Interpretation |
|--------|------------------|----------------|
| Identity | Exact matches | Higher = more similar |
| Conservation | Most common residue frequency | Higher = less variable |
| Entropy | Variability | Lower = more conserved |
| Information Content | Constraint level | Higher = more constrained |

## Tips

- Identity can be calculated different ways (with/without gaps) - specify which you need
- Conservation and entropy are inversely related
- For proteins, use BLOSUM62 for scoring; for DNA, use simple match/mismatch
- Gap-rich columns often indicate alignment uncertainty
- Compare metrics across columns to find functional regions

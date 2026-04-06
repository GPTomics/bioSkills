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

- Percent identity has **four common definitions** (different denominators) producing up to 11.5% difference on the same alignment. Always specify which method: PID1 (gap-inclusive), PID2 (aligned pairs only), PID3 (shorter sequence), PID4 (mean length, recommended for evolutionary studies)
- Conservation and entropy are inversely related; both should be computed ignoring gap characters for interpretable results
- For proteins, use BLOSUM62 for scoring; for DNA, use simple match/mismatch
- Gap-rich columns often indicate alignment uncertainty or guide tree artifacts rather than true biology
- For critical analyses (phylogenetics, selection), quantify alignment confidence per column with GUIDANCE2 or MUSCLE5 ensemble before inference
- Alignment uncertainty propagates to downstream results: different aligners can support different tree topologies; always report alignment method and consider sensitivity analysis
- Average pairwise identity <25% (protein) signals the twilight zone where alignment reliability is questionable and structural methods should be considered

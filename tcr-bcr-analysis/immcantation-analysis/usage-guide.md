# Immcantation Analysis - Usage Guide

## Overview

Analyze BCR repertoires for somatic hypermutation, clonal evolution, and selection using the Immcantation framework (alakazam, shazam, tigger, dowser).

## Prerequisites

```r
install.packages(c('alakazam', 'shazam', 'tigger', 'dowser', 'scoper'))

# For lineage tree building
# Install PHYLIP: http://evolution.genetics.washington.edu/phylip.html
```

## Quick Start

Tell your AI agent:
- "Analyze somatic hypermutation in my BCR data"
- "Build clonal lineage trees"
- "Test for selection in my B cell clones"
- "Infer novel germline alleles"

## Example Prompts

### Mutation Analysis

> "Calculate mutation frequencies for my BCR sequences"

> "Compare SHM levels between naive and memory B cells"

> "What is the R/S ratio in the CDR regions?"

### Clonal Analysis

> "Cluster sequences into clonal lineages"

> "Build phylogenetic trees for expanded clones"

> "Identify the unmutated common ancestor"

### Selection

> "Test for positive selection in my BCR clones"

> "Is there evidence of antigen-driven selection?"

> "Compare selection between conditions"

## What the Agent Will Do

1. Load AIRR-formatted BCR data
2. Cluster sequences into clones (hierarchicalClones)
3. Calculate mutation frequencies (observedMutations)
4. Optionally test selection (BASELINe)
5. Build lineage trees for multi-sequence clones
6. Generate visualizations

## Tips

- **AIRR format** is required - convert from MiXCR if needed
- **Clone threshold 0.15** is typical for IgG; adjust for IgM (0.10)
- **Multiple sequences per clone** needed for trees
- **Selection analysis** requires sufficient clone sizes
- **Germline inference** improves accuracy for individuals

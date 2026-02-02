# Immunogenicity Scoring - Usage Guide

## Overview

Score and rank epitopes for immunogenicity using multi-factor models to prioritize vaccine candidates.

## Prerequisites

```bash
pip install mhcflurry pandas numpy
mhcflurry-downloads fetch
```

## Quick Start

Tell your AI agent what you want to do:
- "Rank these neoantigens by immunogenicity"
- "Score my epitope candidates for vaccine design"
- "Prioritize peptides considering binding and expression"

## Example Prompts

### Scoring

> "Calculate immunogenicity scores for these peptides"

> "Which neoantigens are most likely to be immunogenic?"

### Prioritization

> "Select top 20 vaccine candidates from my list"

> "Rank epitopes considering all immunogenicity factors"

### Analysis

> "Check if these peptides are self-like"

> "Assess anchor residue preferences for HLA-A*02:01"

## What the Agent Will Do

1. Calculate MHC binding scores
2. Assess proteasomal processing likelihood
3. Check self-similarity to proteome
4. Consider expression and clonality
5. Compute weighted immunogenicity score
6. Rank and tier candidates

## Tips

- **Binding weight** - Usually highest (25-35% of score)
- **Agretopicity** - Important for neoantigens (MT vs WT)
- **Self-similarity** - High similarity suggests tolerance
- **Expression** - Unexpressed antigens won't be presented
- **Multiple factors** - No single factor determines immunogenicity

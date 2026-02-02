# Metabolic Reconstruction - Usage Guide

## Overview

Build genome-scale metabolic models from genome sequences using automated reconstruction tools.

## Prerequisites

```bash
# CarveMe (recommended)
pip install carveme
# Requires DIAMOND or BLAST for sequence alignment

# gapseq (alternative)
git clone https://github.com/jotech/gapseq
```

## Quick Start

Tell your AI agent what you want to do:
- "Build a metabolic model from this genome sequence"
- "Reconstruct a model for my bacterial isolate"
- "Create models for all genomes in this directory"

## Example Prompts

### Single Genome

> "Reconstruct a metabolic model from genome.faa"

> "Build a model for this E. coli strain"

### Gap-Filling

> "Create a model that can grow on M9 minimal media"

> "Build and gap-fill a model for LB medium"

### Batch Processing

> "Reconstruct models for all FASTAs in the genomes folder"

> "Build models for my MAGs from metagenome assembly"

## What the Agent Will Do

1. Accept protein FASTA file (or genome for gapseq)
2. Run CarveMe/gapseq reconstruction
3. Gap-fill for specified media if requested
4. Load and inspect the draft model
5. Report basic statistics (reactions, metabolites, genes)
6. Test if model can grow

## Tips

- **Input format** - CarveMe requires protein FASTA (.faa); gapseq can use nucleotide
- **Gap-filling** - Always gap-fill for a defined media to ensure growth
- **Gram stain** - Specify --grampos for Gram-positive bacteria
- **Model size** - Typical bacteria: 1000-2500 reactions
- **Quality check** - Always test growth after reconstruction

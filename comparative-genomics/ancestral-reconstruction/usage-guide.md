# Ancestral Reconstruction - Usage Guide

## Overview

Reconstruct ancestral sequences at phylogenetic nodes using PAML and IQ-TREE for protein resurrection studies and evolutionary trajectory analysis.

## Prerequisites

```bash
# PAML
conda install -c bioconda paml

# IQ-TREE2 (modern alternative)
conda install -c bioconda iqtree

# Python dependencies
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:
- "Reconstruct the ancestral sequence at the root of this tree"
- "Find ancestral states for protein resurrection"
- "Identify ambiguous ancestral positions in my alignment"

## Example Prompts

### Basic ASR

> "Reconstruct ancestral sequences at all internal nodes"

> "Get the ancestral sequence of the last common ancestor"

### Confidence Analysis

> "Which ancestral positions have low confidence?"

> "Show me alternative states at ambiguous positions"

### Protein Resurrection

> "Design constructs for ancestral protein resurrection"

> "Compare ancestral sequence to extant references"

## What the Agent Will Do

1. Prepare alignment in appropriate format
2. Run PAML codeml/baseml or IQ-TREE with ancestral reconstruction
3. Parse RST/state file for ancestral sequences
4. Extract posterior probabilities for each site
5. Identify ambiguous positions with alternative states
6. Calculate overall sequence confidence
7. Suggest alternative constructs for experimental validation

## Tips

- **Confidence threshold** - P > 0.95 is high confidence; P < 0.80 suggests ambiguity
- **Ambiguous sites** - Test alternative states experimentally for resurrection studies
- **Model selection** - Use appropriate substitution model (WAG/LG for proteins)
- **Tree quality** - Ancestral reconstruction depends on tree accuracy
- **Gap handling** - IQ-TREE handles gaps better than PAML
- **Branch length** - Long branches reduce confidence; add intermediate taxa if possible

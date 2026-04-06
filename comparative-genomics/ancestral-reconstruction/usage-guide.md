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

1. Verify tree is rooted (outgroup or midpoint rooting)
2. Select substitution model (LG/WAG+G for proteins; run ModelFinder if unsure)
3. Run PAML codeml/baseml or IQ-TREE with marginal ancestral reconstruction
4. Parse RST/state file for ancestral sequences and per-site probabilities
5. Assess overall reconstruction quality (fraction of high-confidence sites)
6. Identify ambiguous positions with plausible alternative states
7. Flag sites where epistatic interactions may affect function
8. Suggest primary ML construct plus alternatives at ambiguous positions

## Tips

- **Rooted tree required** - ASR depends on root position; use outgroup or midpoint rooting
- **Joint vs marginal** - Use marginal for resurrection studies (gives per-site probabilities); joint is faster but less informative
- **Confidence threshold** - P > 0.95 is high confidence; P < 0.80 suggests ambiguity; always report ambiguous sites
- **Model selection** - Use LG or WAG with +G for proteins; run ModelFinder for data-driven choice; use same model for tree and ASR
- **Epistasis** - ML ancestral sequence may be non-functional due to untested residue combinations; test multiple alternative constructs
- **Alignment** - Use PRANK for coding sequences (correctly handles insertions vs deletions)
- **Gap handling** - Gaps are indels, not substitutions; PAML `cleandata=0` treats as ambiguity, `cleandata=1` removes gapped columns; IQ-TREE handles gaps better
- **Taxon sampling** - Dense sampling near the node of interest improves accuracy; add intermediate taxa to break long branches
- **Depth limits** - ASR becomes unreliable for very deep divergences; proteins are better than nucleotides for ancient nodes

## Related Skills

- comparative-genomics/positive-selection - Selection analysis on ancestral branches
- comparative-genomics/ortholog-inference - Identify orthologs for reconstruction
- phylogenetics/modern-tree-inference - Generate rooted trees for ASR
- alignment/multiple-alignment - PRANK alignment (indel-aware)

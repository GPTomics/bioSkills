# Positive Selection - Usage Guide

## Overview

Detect positive selection using dN/dS tests with PAML codeml and HyPhy to identify sites and branches under adaptive evolution.

## Prerequisites

```bash
# PAML
conda install -c bioconda paml

# HyPhy (recommended for episodic selection)
conda install -c bioconda hyphy

# Python dependencies
pip install biopython scipy
```

## Quick Start

Tell your AI agent what you want to do:
- "Test for positive selection on this gene across mammals"
- "Find positively selected sites in my protein alignment"
- "Run a branch-site test for selection on the primate lineage"

## Example Prompts

### Site-Level Selection

> "Find sites under positive selection in this immune gene"

> "Run PAML M8 vs M7 test on my codon alignment"

### Branch-Specific Selection

> "Test if the human branch shows positive selection"

> "Run branch-site test marking primates as foreground"

### HyPhy Analysis

> "Use BUSTED to test for gene-wide positive selection"

> "Run MEME to find episodically selected sites"

## What the Agent Will Do

1. Prepare codon alignment in PHYLIP format
2. Create codeml control file for appropriate model
3. Run site models (M7 vs M8) or branch-site tests
4. Perform likelihood ratio test for significance
5. Extract positively selected sites (BEB posterior > 0.95)
6. Alternatively run HyPhy BUSTED/MEME for episodic selection
7. Report sites with significance levels

## Tips

- **Codon alignment** - Must be in-frame; remove sequences with frameshifts
- **Site models** - Use M8 vs M7 for most analyses; M8 vs M8a is more stringent
- **Branch-site** - Mark foreground with #1 in tree; useful for lineage-specific selection
- **BEB thresholds** - P > 0.95 significant (*), P > 0.99 highly significant (**)
- **HyPhy preference** - Use MEME for episodic selection, FEL for pervasive
- **Multiple testing** - Correct p-values when testing many genes

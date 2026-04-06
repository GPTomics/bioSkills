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

1. Align CDS sequences with PRANK (codon-aware) or verify existing alignment quality
2. Screen for recombination with GARD; partition alignment if breakpoints found
3. Run gene-wide selection screen (BUSTED) to test for any positive selection
4. If significant, run site models (M8 vs M8a or M7 vs M8) and/or MEME for site identification
5. Perform likelihood ratio test with appropriate chi-squared distribution
6. Extract positively selected sites (BEB posterior > 0.95)
7. Check for false-positive signals (saturation, gBGC, alignment artifacts at flagged sites)
8. Apply multiple testing correction (FDR) when analyzing gene families

## Tips

- **Alignment first** - Use PRANK for codon alignment (correctly handles insertions); avoid Gblocks filtering
- **Check for recombination** - Run GARD before any selection test to avoid inflated false positives
- **Site models** - M8 vs M8a is the recommended primary test (more stringent than M7 vs M8)
- **Branch-site LRT** - Uses 50:50 chi-squared mixture distribution; critical value is 2.71 (not 3.84)
- **BEB thresholds** - P > 0.95 significant (*), P > 0.99 highly significant (**)
- **HyPhy pipeline** - BUSTED (gene-wide) -> aBSREL (branch) -> MEME (episodic sites) / FEL (pervasive sites)
- **False positives** - dN/dS > 1 can result from GC-biased gene conversion, saturation (dS > 3), or alignment errors rather than positive selection
- **Multiple testing** - Use FDR (Benjamini-Hochberg) across genes, not Bonferroni (genes are non-independent)
- **Saturation check** - If dS > 3, codon-based analysis is unreliable; switch to amino acid comparison

## Related Skills

- comparative-genomics/synteny-analysis - Synteny context for gene pairs
- comparative-genomics/ortholog-inference - Identify orthologs for analysis
- comparative-genomics/ancestral-reconstruction - Reconstruct sequences at selected branches
- alignment/msa-parsing - Parse and manipulate codon alignments
- alignment/multiple-alignment - PRANK codon-aware alignment
- phylogenetics/modern-tree-inference - Generate trees for codeml

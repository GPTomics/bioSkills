# Species Trees - Usage Guide

## Overview

Estimate species-level phylogenies from multi-locus data using coalescent methods that account for gene tree discordance caused by incomplete lineage sorting (ILS). Covers ASTRAL-III (and the broader ASTER suite), SVDQuartets, BPP, concordance factor interpretation, the anomaly zone, and deciding when concatenation is appropriate versus when coalescent methods are required.

## Prerequisites

```bash
# ASTER package (includes ASTRAL-III, wASTRAL, ASTRAL-Pro, CASTER)
# Build from source: https://github.com/chaoszhang/ASTER
git clone https://github.com/chaoszhang/ASTER.git
cd ASTER
make

# IQ-TREE2 (for gene trees and concordance factors)
conda install -c bioconda iqtree

# PAUP* (for SVDQuartets)
# Download from https://paup.phylosolutions.com/

# BPP (for species delimitation)
conda install -c bioconda bpp
```

## Quick Start

Tell your AI agent what you want to do:
- "Infer a species tree from my multi-locus dataset using ASTRAL"
- "Compute concordance factors to assess gene tree agreement with my species tree"
- "Should I concatenate or use coalescent methods for my phylogenomic dataset?"
- "Run SVDQuartets on my RADseq data since per-locus alignments are too short for gene trees"
- "Check whether gene tree discordance in my dataset is from ILS or introgression"

## Example Prompts

### Species Tree Estimation
> "I have 200 locus alignments in loci/. Infer gene trees and then estimate a species tree with ASTRAL"

> "Run wASTRAL instead of ASTRAL-III on my gene trees to account for gene tree estimation uncertainty"

> "Estimate a species tree from my multi-copy gene family trees using ASTRAL-Pro"

### Concordance Factor Analysis
> "Compute gene concordance factors and site concordance factors for my species tree"

> "My species tree has gCF below 30% at several backbone nodes. What does that mean?"

> "Compare gCF and sCF values to determine whether discordance is biological or from estimation error"

### Concatenation vs Coalescent Decision
> "My phylogenomic dataset has 500 loci. Should I concatenate or use coalescent methods?"

> "I ran both concatenation and ASTRAL and got different topologies. Which should I trust?"

> "Are my loci in the anomaly zone? How do I check?"

### SVDQuartets and Alternative Methods
> "Run SVDQuartets on my UCE dataset in PAUP*"

> "My per-locus alignments are only 300 bp. What species tree method should I use?"

> "Run CASTER directly on my alignments without inferring gene trees first"

### Species Delimitation
> "Use BPP to jointly estimate the species tree and test species boundaries for my populations"

> "Run BPP A01 analysis to test species delimitation on my guide tree"

### Introgression Detection
> "Check quartet asymmetry in my gene trees to test for introgression"

> "Gene tree discordance is too high to be explained by ILS alone. What network methods should I try?"

## What the Agent Will Do

1. Assess the dataset: number of loci, alignment lengths, taxon sampling, and whether gene trees already exist
2. Infer per-locus gene trees with IQ-TREE2 if not already available
3. Compute concordance factors to quantify discordance levels
4. Evaluate whether concatenation is appropriate based on gCF distribution
5. Run ASTRAL (or wASTRAL for noisy gene trees) for species tree estimation
6. Interpret branch support (local posterior probabilities) and concordance factors
7. Check quartet asymmetry for introgression signal
8. Recommend SVDQuartets or CASTER when per-locus alignments are too short
9. Suggest BPP for species delimitation with small taxon sets

## Tips

- Always compute concordance factors before deciding between concatenation and coalescent methods. The gCF distribution is the key diagnostic
- ASTRAL branch lengths are in coalescent units, not substitutions per site; short branches (<0.5 CU) predict high ILS
- ASTRAL uses local posterior probabilities, not bootstrap; the threshold for strong support is >= 0.95
- wASTRAL outperforms ASTRAL-III when gene trees are noisy (short loci, low bootstrap); it uses branch length and support information
- For RADseq or UCE data with short per-locus alignments, SVDQuartets or CASTER avoids the gene tree estimation bottleneck
- The anomaly zone is not rare; rapid radiations commonly fall within it, making concatenation statistically inconsistent
- If one minor quartet is much more frequent than the other (q2 >> q3), suspect introgression rather than ILS
- BPP is the gold standard for species delimitation but is computationally expensive; reserve it for small datasets (<20 species, <50 loci)
- The ASTER C++ package replaces the original Java ASTRAL implementations with significantly better performance

## Related Skills

- modern-tree-inference - ML gene tree inference needed as ASTRAL input
- bayesian-inference - Bayesian tree inference; StarBEAST3 for species trees with dating
- phylogenetics/divergence-dating - Dating divergences after species tree estimation
- tree-manipulation - Rooting and modifying species trees
- tree-visualization - Visualizing species trees with concordance annotations
- alignment/multiple-alignment - Alignment quality affects gene tree accuracy

# Species Trees - Usage Guide

## Overview

This skill estimates a species tree from multi-locus phylogenomic data under the multispecies coalescent (MSC). The central idea is that a species tree is not a gene tree: each locus has its own genealogy, and under the MSC those genealogies disagree with the species tree (and with each other) by incomplete lineage sorting (ILS) alone, even with zero estimation error. That discordance is signal, not noise. The consequence that breaks intuition is the anomaly zone (two short successive speciation events): the single most common gene tree can be the wrong species tree, and concatenation (which pretends all loci share one tree) is statistically inconsistent and positively misleading there, converging on the wrong topology with full support as more loci are added. The support is the trap.

Coalescent summary methods (the ASTRAL family) stay consistent because they ask which species tree is most consistent with the observed gene-tree distribution, working at the quartet level where no anomaly zone exists. The deliverable is concordance, not a single best method: run both concatenation and a coalescent method, compute gene and site concordance factors, and trust the coalescent tree exactly where the two disagree on low-gCF branches.

## Prerequisites

```bash
# ASTER package (provides astral/ASTRAL-III, wastral, astral-pro)
conda install -c bioconda aster      # or build from https://github.com/chaoszhang/ASTER

# IQ-TREE2 (per-locus gene trees and concordance factors)
conda install -c bioconda iqtree

# Optional, for short-locus and small-dataset methods
# PAUP* (SVDQuartets): https://paup.phylosolutions.com/
conda install -c bioconda bpp        # BPP (delimitation, full-likelihood)
```

- Conceptual: a species tree is a distribution-of-genealogies problem, not one tree problem. Know your approximate ILS level (internode lengths in coalescent units, ancestral Ne) before choosing a method.
- Input: per-locus gene trees (one Newick per line) for the ASTRAL family, or alignments for SVDQuartets/BPP. Gene trees should carry branch support so weak branches can be contracted or weighted.
- Orthology: standard ASTRAL needs single-copy gene trees; multi-copy families need ASTRAL-Pro (see comparative-genomics/ortholog-inference).

## Quick Start

Tell your AI agent what you want to do:
- "Infer a species tree from my multi-locus gene trees with ASTRAL"
- "Should I concatenate or use a coalescent method for my phylogenomic dataset?"
- "Compute gene and site concordance factors for my species tree and tell me which branches are contested"
- "My per-locus alignments are only 300 bp: which species-tree method avoids the gene-tree bottleneck?"
- "Is the discordance in my data ILS or introgression?"

## Example Prompts

### Species Tree Estimation
> "I have 200 locus alignments in loci/. Infer gene trees and estimate a species tree, contracting weak gene-tree branches first"

> "Run wASTRAL instead of plain ASTRAL because my gene trees are noisy, and explain why"

> "Estimate a species tree from my multi-copy gene-family trees using ASTRAL-Pro without pre-orthology"

### Concatenation vs Coalescent Decision
> "My 500-locus dataset spans a rapid radiation: is concatenation safe or am I in the anomaly zone?"

> "I ran concatenation and ASTRAL and got different backbone topologies: which should I trust and how do I decide?"

### Concordance and Discordance
> "Compute gCF and sCF for my species tree and flag branches with high support but low concordance"

> "My backbone nodes have gCF below 30: what does that mean and is the tree resolved there?"

> "Check the q2/q3 quartet symmetry on my contested branches to tell ILS apart from introgression"

### Short Loci and Small Datasets
> "Run SVDQuartets on my UCE dataset because per-locus gene trees are unreliable"

> "Use BPP to jointly estimate the species tree and test species boundaries for my populations"

## What the Agent Will Do

1. Assess the dataset: number of loci, alignment lengths, taxon sampling, approximate ILS level, and whether gene trees already exist.
2. Infer per-locus gene trees with support if needed (routes to modern-tree-inference), then contract branches below ~10% support to polytomies.
3. Run wASTRAL as the primary estimate, or `astral` for the classic localPP / polytomy-test workflow, mapping the correct ASTER vs Java flags (`-u` annotation in ASTER, `-t` annotation in Java).
4. Run concatenation in parallel as a contrast, not as the answer.
5. Compute gene and site concordance factors against the species tree and map discordance onto branches.
6. Read localPP on its own scale (not as bootstrap) and use the polytomy test on weak branches.
7. Inspect q2/q3 symmetry on contested branches; if asymmetric, route to introgression methods rather than forcing a tree.
8. Recommend SVDQuartets (site-based) for short loci, or BPP/StarBEAST2 for small high-stakes datasets, delimitation, or dating.

## Tips

- A species tree is not a gene tree: treat the data as a distribution of genealogies, and expect discordance even with perfect data.
- More loci do not rescue concatenation in the anomaly zone; they make the wrong tree more confident. Measure ILS first.
- ASTRAL branch lengths are in coalescent units, not time or substitutions, and tip lengths are undefined; for dates use StarBEAST2 or divergence-dating.
- localPP is a coalescent posterior, not a bootstrap. localPP = 1.0 is expected on resolved branches; localPP ~ 0.33 is a three-way tie.
- Report gCF/sCF on every phylogenomic tree; bootstrap and localPP saturate at scale and hide the conflict that concordance factors expose.
- Contract gene-tree branches below ~10% support before ASTRAL, or use wASTRAL, which weights continuously and generally beats hard contraction.
- Symmetric minority quartets (q2 ~ q3) mean ILS; asymmetric (q2 != q3, significant D / HyDe / QuIBL) mean introgression; confirm the asymmetry before invoking gene flow.
- In ASTER `astral`, `-t` is threads and `-u` is annotation; in classic Java ASTRAL `-t` is the annotation level. Confirm which is installed before scripting.

## Related Skills

- modern-tree-inference - per-locus ML gene trees (ASTRAL input) and gCF/sCF computation
- bayesian-inference - full Bayesian co-estimation; StarBEAST2 for species tree plus dates
- divergence-dating - turning a coalescent-unit species tree into a dated tree
- tree-io - reading and writing the Newick gene-tree files ASTRAL consumes
- comparative-genomics/ortholog-inference - single-copy vs multi-copy decision upstream (ASTRAL vs ASTRAL-Pro)

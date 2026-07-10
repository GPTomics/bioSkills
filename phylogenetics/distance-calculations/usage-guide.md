# Distance Calculations - Usage Guide

## Overview

This skill computes evolutionary distance matrices from sequence alignments and builds fast distance-based trees (NJ, BIONJ, FastME, UPGMA), with bootstrap support and a saturation pre-flight. The idea everything hinges on: a distance is a model-corrected estimate of expected substitutions per site, the correction is where the biology lives, and the matrix discards the per-site information ML uses, so distance trees are fast and scalable but only as good as the correction, and structurally weaker than ML on hard, homoplasy-rich data. Two flat rules drive most decisions: never use UPGMA for molecular phylogeny without an established clock, and never call Biopython's `identity` distance a Jukes-Cantor tree (it is uncorrected; use ape `dist.dna` or FastME for a real correction).

## Prerequisites

```bash
pip install biopython scikit-bio        # Python NJ algorithm + p-distance / matrix distances
# Model-corrected distances live in R:
# Rscript -e "install.packages(c('ape','phangorn'))"
# FastME standalone (best distance tree, large n): conda install -c bioconda fastme
```

Conceptual prerequisites: a trimmed multiple sequence alignment (distances are only as good as the alignment); awareness that Bio.Phylo and scikit-bio do NOT correct for multiple hits, so divergent data needs ape/FastME for the distance step.

## Quick Start

Tell your AI agent what you want to do:
- "Build a neighbor joining tree from this alignment"
- "Compute a TN93-corrected distance matrix and a FastME tree"
- "Test my sequences for substitution saturation before I trust a distance tree"
- "Make a bootstrap consensus distance tree with 500 replicates"

## Example Prompts

### Distance Matrices
> "Compute a model-corrected distance matrix from this DNA alignment using TN93 with a gamma correction"

> "My taxa have very different GC content: give me a LogDet distance matrix instead of JC"

> "Calculate protein distances from this alignment using the LG matrix"

### Tree Building
> "Build a FastME balanced-minimum-evolution tree from these sequences"

> "Make a neighbor joining tree and explain why you did not use UPGMA"

> "Generate a fast NJ starting tree to seed an IQ-TREE ML search"

### Saturation and Diagnostics
> "Check whether these deep divergences are saturated before building a distance tree"

> "Plot transitions against corrected distance to see if the signal has plateaued"

### Bootstrap Support
> "Build a distance tree and report bootstrap support with 1000 replicates"

> "Why is my bootstrap support high on a clade that looks wrong?"

## What the Agent Will Do

1. Read and sanity-check the alignment (flag ambiguous blocks to trim upstream)
2. Choose a correction from the data: shallow/barcoding -> p-distance or K80; divergent DNA -> TN93 (+ gamma if ASRV); compositional skew -> LogDet; protein -> LG/WAG
3. Compute the distance matrix in the right tool (ape/FastME for real corrections; Bio.Phylo/scikit-bio only for identity/score distances)
4. Run a saturation pre-flight (Xia Iss vs Iss.c) for deep data and refuse to over-trust a saturated tree
5. Build the tree (FastME for the best distance tree; NJ/BIONJ for speed or ML seeding; UPGMA only if a clock is established)
6. Bootstrap if support is wanted, and report it as precision, not accuracy

## Tips

- Distance buys speed and scale by discarding information; that trade is fine on shallow/clean data or huge n, and a mistake on deep, homoplasy-rich data where the discarded signal is the whole point; route those to modern-tree-inference.
- NJ is statistically consistent only if the distances are correct; a consistent algorithm fed saturated distances converges confidently on the wrong tree.
- If a topology flips between JC and LogDet, suspect compositional heterogeneity, not noise.
- Estimating gamma alpha usually needs a preliminary ML tree, so a distance pipeline that needs ASRV is already leaning on ML.
- Biopython's `DistanceCalculator` has no JC/K80/TN93; do the correction in ape/FastME and pass the matrix in, or stay in R.
- 100-1000 bootstrap replicates is standard; clades below ~70% are conventionally unsupported.

## Related Skills

- modern-tree-inference - ML inference and where a distance tree seeds the ML search
- tree-manipulation - rooting and pruning the unrooted trees these methods emit
- tree-io - saving and converting constructed trees
- alignment/alignment-io - reading and trimming the alignment that gates every distance

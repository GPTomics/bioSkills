# Immcantation Analysis - Usage Guide

## Overview

Immcantation analyzes B-cell receptor repertoires: it partitions somatically-hypermutated sequences into clonal families, quantifies somatic hypermutation (SHM) and antigen-driven selection, and reconstructs antibody lineage trees. The step everything inherits from is the clonal-clustering threshold: `shazam::distToNearest` builds a bimodal distance-to-nearest-neighbor distribution and `shazam::findThreshold` locates its valley, which becomes the per-dataset cutoff that every downstream number (clone counts, diversity, selection, trees) inherits. Because B cells hypermutate, exact-CDR3 clonotypes are wrong for BCR; clones are grouped by shared V gene, J gene, and junction length, then clustered by junction nucleotide distance. Germline reconstruction and TIGGER genotyping must precede any mutation counting, and diversity must be compared at equal sampling depth.

## Prerequisites

```r
install.packages(c('alakazam', 'shazam', 'tigger', 'dowser', 'scoper'))
# Recommended: run inside the immcantation/suite Docker image, which bundles
# IgBLAST, Change-O, IMGT germline setup, PHYLIP, and IgPhyML (the fragile installs).
```

External dependencies (for upstream annotation and trees): IgBLAST + Change-O (`AssignGenes.py`, `MakeDb.py`) to produce the AIRR TSV, IMGT-gapped germline references, and IgPhyML/PHYLIP for lineage trees.

## Quick Start

Tell your AI agent what you want to do:
- "Derive the clonal threshold for my BCR data and assign clones"
- "Quantify somatic hypermutation by CDR and FWR region"
- "Test for antigen-driven selection in my B cell clones"
- "Personalize the germline genotype with TIGGER first"
- "Build IgPhyML lineage trees for my expanded clones"

## Example Prompts

### Clonal Assignment

> "Compute distToNearest, find the threshold from the bimodal histogram, and cluster my IGH sequences into clones"

> "My distance histogram is unimodal - cluster with spectralClones instead"

> "This is single-cell BCR - cluster on heavy chains and split by light chain"

### Mutation and Selection

> "Reconstruct the D-masked germline and measure R and S mutation frequency in CDR and FWR"

> "Run BASELINe and tell me whether the CDRs show positive selection"

> "Genotype the subject with TIGGER before counting mutations so polymorphisms are not scored as SHM"

### Diversity and Lineages

> "Compare clonal diversity between timepoints at equal sequencing depth"

> "Build lineage trees for clones with at least three sequences and color tips by isotype"

> "Reconstruct the ancestral antibody at the root of the largest clone"

## What the Agent Will Do

1. Load the AIRR-formatted BCR table and confirm required columns are present.
2. Personalize the germline genotype with TIGGER (novel alleles, `inferGenotypeBayesian`, `reassignAlleles`).
3. Reconstruct the per-sequence germline with `createGermlines`.
4. Compute `distToNearest`, derive the threshold with `findThreshold`, and inspect the histogram.
5. Cluster sequences into clones with `hierarchicalClones` (or `spectralClones` if unimodal).
6. Rebuild the per-clone germline and measure SHM with `observedMutations` (V region only).
7. Optionally test selection with `calcBaseline`/`groupBaseline`.
8. Compare Hill-number diversity at equal depth with `alphaDiversity`.
9. Build IgPhyML lineage trees with `formatClones`/`getTrees` and plot them.

## Tips

- The threshold is per-dataset: never reuse a number from a paper. Re-derive it for each study, locus, and depth.
- Run TIGGER first: an unrecorded personal allele reads as recurrent SHM and also corrupts `distToNearest`.
- Mask the junction: use the D-masked germline and `regionDefinition = IMGT_V`; junctional bases have no germline template.
- S5F is a targeting model (`HH_S5F`), used for selection nulls, not a `mutationDefinition`. There is no `MUTATION_SCHEMES$S5F`.
- Report mutation frequency, not raw counts, when sequence coverage or length varies.
- Single-cell: cluster on the heavy chain, then split clones by light chain with `dowser::resolveLightChains` (the old scoper `only_heavy`/`split_light` args are deprecated); light chains alone cannot define clones.
- Never pool clones across individuals; genotypes and novel alleles are private.
- Compare diversity only after uniform resampling to equal N (the `alphaDiversity` default) with bootstrap CIs.
- Prefer Dowser/IgPhyML over legacy `buildPhylipLineage`: standard phylogenetic models violate SHM context-dependence, non-reversibility, and the known germline root.

## Related Skills

- mixcr-analysis - Produce AIRR/clonotype input for BCR
- scirpy-analysis - Single-cell BCR integration and handoff
- specificity-annotation - Convergent/public antibody signatures
- phylogenetics/tree-visualization - General lineage-tree plotting concepts
- phylogenetics/modern-tree-inference - Phylogenetic inference background
- workflows/tcr-pipeline - End-to-end orchestration

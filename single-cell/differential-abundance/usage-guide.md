# Differential Abundance Testing - Usage Guide

## Overview

Differential abundance testing asks whether cell-type proportions or composition changed between conditions (treatment vs control, disease vs healthy). Because proportions live on a simplex (they sum to 1 and are not independent), naive per-cluster proportion tests are invalid. This skill covers cluster-free neighborhood testing (Milo) and cluster-based compositional models (scCODA, sccomp, propeller), and how to keep compositional shifts from being misread as differential expression.

## Prerequisites

```r
# Milo, sccomp, propeller (R)
BiocManager::install(c('miloR', 'sccomp', 'speckle'))
```

```bash
# scCODA (Python)
pip install scanpy sccoda
```

## Quick Start

Tell your AI agent what you want to do:
- "Did any cell-type proportions change between my conditions?"
- "Test differential abundance with Milo on my integrated data"
- "Run scCODA to find populations that expanded with treatment"
- "Check whether my DE signal is actually a compositional shift"

## Example Prompts

### Cluster-free abundance
> "Run Milo differential abundance on my kNN graph and report SpatialFDR neighborhoods"
> "Find transitional states that expanded with treatment without committing to clusters"

### Cluster-based composition
> "Test cell-type proportion changes with scCODA using a stable reference cell type"
> "Run sccomp on my count table, robust to outlier samples"
> "Use propeller to quickly test proportion differences across groups"

### Guarding the DE/DA confound
> "Pair my pseudobulk DE with a differential-abundance test"
> "Is the change in this cluster's expression real or just a shift in substate proportions?"

### Interpretation
> "Which population actually drives the compositional change relative to the reference?"
> "Do Milo and scCODA agree on which cell types changed?"

## What the Agent Will Do

1. Confirm there are biological replicates per condition (samples, not cells, are the unit)
2. Choose cluster-free (Milo) and/or cluster-based (scCODA/sccomp/propeller) testing for the question
3. Build a per-sample cell-type count table or a Milo object from an integrated embedding
4. For scCODA, select a stable reference cell type (or automatic) and explain the relative interpretation
5. Fit the model and report FDR/SpatialFDR-controlled effects with direction and uncertainty
6. Reconcile cluster-free and cluster-based results when both are run
7. Cross-check against the condition-DE analysis to separate abundance shifts from expression changes

## Tips

- **Samples are the replicate** - more donors increase power for an abundance claim; more cells per donor barely do.
- **Never run a per-cluster proportion t-test** - the simplex constraint makes one expansion force spurious depletions elsewhere.
- **Choose the scCODA reference deliberately** - a reference that actually changes biases every other call; use automatic if unsure.
- **Report SpatialFDR for Milo** - it corrects for overlapping neighborhoods; raw p-values over-call.
- **Milo k and prop trade resolution for power** - larger neighborhoods are better powered but blur fine shifts.
- **Always pair DA with DE** - a compositional shift can masquerade as differential expression and vice versa.
- **Build Milo's graph on an integrated embedding** - batch structure in the kNN graph produces false abundance calls.

## Related Skills

- clustering - Define the clusters whose abundance is tested
- cell-annotation - Annotate cell types before testing their proportions
- markers-annotation - Pair condition DE with abundance testing to separate the confound
- batch-integration - Build the integrated embedding Milo's kNN graph relies on
- differential-expression/deseq2-basics - Pseudobulk condition DE that abundance testing complements
- pathway-analysis/go-enrichment - Characterize the populations that expanded or contracted

# Spatial Cell-Cell Communication - Usage Guide

## Overview

This skill infers cell-cell communication from spatial transcriptomics data by scoring ligand-receptor co-expression between sender and receiver cell types against a permutation null. Its core job is decision-making, not recipe-running: choosing a method by whether spatial distance is actually modeled (squidpy ligrec is space-blind, COMMOT is distance-aware, stLearn and CellChat-spatial sit between), choosing the ligand-receptor database deliberately because it drives the result as much as the algorithm, guarding against segmentation spillover that fabricates short-range hits, and treating every call as a co-expression hypothesis on a confidence ladder rather than as validated signaling.

## Prerequisites

```bash
pip install squidpy scanpy anndata statsmodels matplotlib
pip install commot          # distance-aware optimal-transport inference
# CellChat v2 (spatial) and NicheNet are R/Bioconductor packages
```

## Quick Start

Tell your AI agent what you want to do:
- "Map cell-cell communication in my Visium data and tell me which calls are spatially supported"
- "Run ligand-receptor analysis but treat the hits as hypotheses, not signaling"
- "Use a distance-aware method so adjacency is actually modeled"
- "Check whether a short-range interaction is just segmentation spillover"

## Example Prompts

### Method Choice
> "I have Xenium data -- pick a communication method that actually models spatial distance and explain why squidpy ligrec alone is not spatial."

> "Run COMMOT on my spatial data and set the signaling range correctly for a secreted chemokine versus a contact-dependent Notch pair."

### Honest Inference
> "Run ligand-receptor analysis and correct for multiple testing across all pair and cell-type-pair combinations."

> "Rank the communication hits as hypotheses and tell me what evidence would move each up the confidence ladder."

### Auditing
> "A T-cell-to-macrophage interaction came up at their shared boundary -- check whether this is transcript spillover from imperfect segmentation."

> "My 300-gene panel returned no communication -- is that a real negative or a panel limitation?"

## What the Agent Will Do

1. Confirm cell-type annotations exist and decide whether the platform even supports the analysis (panel coverage of ligands/receptors).
2. Choose a method by whether spatial distance must be modeled and by secreted-vs-contact range.
3. Build or reuse the spatial neighbor graph in the correct coordinate units.
4. Run the inference (squidpy ligrec baseline and/or distance-aware COMMOT), naming the database and version.
5. Correct for multiple testing over the full pair-by-cell-type-pair space.
6. Audit short-range hits for segmentation spillover and place survivors on the confidence ladder.

## Tips

- A ligand-receptor score is a co-expression HYPOTHESIS, not signaling -- there is no flux, binding, or causality in any standard tool. Refuse "A signals to B" phrasing.
- squidpy `ligrec` is space-blind by default (it permutes cluster labels). Do not present it as spatial communication; use COMMOT or stLearn for distance-modeled inference.
- The L-R database drives the result as much as the algorithm (Dimitrov 2022). Report the database and version; two tools "agreeing" often just share a database.
- Segmentation spillover fabricates short-range co-expression (Mitchel 2026). Validate any boundary-localized hit against segmentation quality before believing it.
- Secreted and contact-dependent ligands need different ranges. One global distance cutoff is wrong for at least one class -- interrogate any fixed radius.
- Correct for the thousands of pair-by-cell-type-pair tests; nominal-p "hotspots" manufacture networks.
- A targeted imaging panel rarely contains the relevant ligands/receptors. A "no communication" result on a small panel is uninformative.
- `threshold` in `sq.gr.ligrec` is the expression-fraction floor, not a p-value cutoff.

## Related Skills

- spatial-transcriptomics/image-analysis - the segmentation-spillover circularity that fabricates short-range hits
- spatial-transcriptomics/spatial-neighbors - build the spatial graph distance-aware methods inherit
- spatial-transcriptomics/spatial-statistics - neighborhood enrichment and co-occurrence for which-types-co-occur questions
- spatial-transcriptomics/spatial-domains - annotate sender and receiver cell types
- single-cell/cell-communication - the non-spatial CellPhoneDB/CellChat/NicheNet baseline
- pathway-analysis/go-enrichment - enrich downstream receiver-response programs

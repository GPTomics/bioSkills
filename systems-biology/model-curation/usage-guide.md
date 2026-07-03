# Model Curation - Usage Guide

## Overview

Curation turns a draft genome-scale model into a trustworthy one, and it splits into two axes that must not be confused. The first is consistency and annotation - mass and charge balance, no stoichiometric leaks, SBO terms, SBML/FBC conformance - which memote scores. The second is predictive validity - does the model reproduce measured growth, carbon-source usage, and gene essentiality - which memote does not score. A model can reach a 90% memote score and still mispredict every knockout, so a high score means well-formed and well-annotated, never "biologically correct." The single most dangerous defect a good score can hide is an energy-generating cycle that makes ATP from nothing; it passes mass balance, inflates growth, and invalidates flux predictions, so test for it explicitly. Gap-filling forces growth on a chosen medium and its added reactions are the least-evidenced part of the model.

## Prerequisites

```bash
pip install memote cobra
# A universal reaction database (e.g. a BiGG universal model SBML) is needed for gap-filling.
```

Inputs: a genome-scale model in SBML (a draft from metabolic-reconstruction, or a published model), and for gap-filling a universal reaction database in the same namespace. For validation, measured growth/essentiality data on a defined medium.

## Quick Start

Tell your AI agent:
- "Run memote on my model and tell me which tests fail, not just the score"
- "Check whether my model can make ATP from nothing"
- "Gap-fill my model to grow on M9, and list what was added so I can flag it"
- "Find unbalanced reactions and dead-end metabolites"
- "My model scores 85% but mispredicts knockouts - what should I actually validate?"

## Example Prompts

### Interpreting a Score
> "Run memote on my model, report the total score, and then explain what the score does and does not tell me about whether the model predicts biology."

### Energy-Generating Cycles
> "Close all exchanges and check whether my model can still produce ATP; if it can, help me find and constrain the reactions forming the erroneous energy-generating cycle."

### Gap-Filling
> "My draft cannot grow on M9 glucose. Gap-fill it from the BiGG universal model, show me the alternative gap-fill sets, and mark the added reactions as low-confidence."

### Validation
> "My model has a high memote score. Compare its predicted gene essentiality to the Keio set on M9 and its carbon-source usage to Biolog data, and report where it fails."

## What the Agent Will Do

1. Load the model and run memote (`memote run` / `report snapshot`), reading which tests fail, not just the total.
2. Run the energy-generating-cycle test (max ATP with all uptake closed must be ~0) and constrain directionality if it fails.
3. Find mass/charge-unbalanced reactions (skip exchange/sink/demand) and fix protons/charge first.
4. Locate dead-end metabolites and decide whether to add a reaction or fix stoichiometry.
5. Gap-fill deliberately toward the objective on a stated medium, recording and flagging added reactions.
6. Validate predictive quality against measured growth/essentiality on the matched medium - the step memote does not perform.

## Tips

- A memote score measures consistency and annotation, not biological correctness; always read which tests fail and validate predictions separately.
- Test for energy-generating cycles explicitly: close all exchanges, maximize the ATP maintenance reaction; a positive result is a defect that inflates growth.
- `gapfill` has no `demand` argument; set the objective and use `lower_bound`/`demand_reactions=False`. Gap-filled reactions are hypotheses, not evidence.
- Skip exchange/sink/demand reactions when checking mass balance - they are intentionally unbalanced.
- Proton (H) and charge imbalance at pH 7 is the most common real balancing fix.
- Dead-end metabolites signal a missing reaction, wrong stoichiometry, or a missing transport/exchange.
- memote scores are only comparable within a memote version; state the version when reporting a score.
- The Python entry points are `memote.suite.api.test_model` and `snapshot_report`, not `run`/`snapshot`.
- Validate against Biolog carbon-source phenotypes and an essentiality screen on the matched medium; report MCC for essentiality (minority class).
- Automation (reconstruction + memote-in-CI + auto-gap-fill) gets you a draft and a QC dashboard, not a validated model; directionality, GPR, biomass, and evidence-based gap-filling stay expert judgments.

## Related Skills

- systems-biology/metabolic-reconstruction - Produces the draft this skill curates
- systems-biology/flux-balance-analysis - Test the curated model's predictions
- systems-biology/gene-essentiality - Validate curation against measured essentiality
- pathway-analysis/kegg-pathways - Source KEGG annotations for reactions/metabolites
- database-access/uniprot-access - Cross-reference gene/protein annotations

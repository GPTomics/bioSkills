# Retrosynthesis Usage Guide

## Overview

Plan synthetic routes from target molecules back to commercial building blocks using maintained AiZynthFinder releases and, when useful, a maintained or exactly version-pinned template-free model. Score routes by configured scorer, reaction count, solved state, building-block availability, and independently configured forward-prediction checks.

## Prerequisites

```bash
pip install aizynthfinder
```

Download AiZynthFinder USPTO model + building-block stocks.

Chemformer is archived legacy research code, not a supported `pip install chemformer` dependency. Reproduction requires the archived repository's documented environment, checkpoint, tokenizer, and configuration at a pinned commit.

## Quick Start

Tell the AI agent what to do:
- "Plan retrosynthesis for this target SMILES using AiZynthFinder"
- "Batch retrosynthesize 100 generated molecules; report feasibility"
- "Find shorter routes using documented AiZynthFinder rewards or post-search route scorers"
- "Validate retrosynthesis with forward prediction via Molecular Transformer"
- "Estimate synthesis cost using building-block prices"

## Example Prompts

### Single target retrosynthesis
> "Plan retrosynthesis for SMILES 'CC(=O)Nc1ccc(C(=O)Nc2ccccc2)cc1' using AiZynthFinder. USPTO templates, ZINC stock. Output top 5 routes with depth, score, building-block list."

### Batch feasibility screening
> "Batch retrosynthesize 100 generated SMILES from candidates.csv. Report whether any route is solved under the configured stock and search budget, plus the shortest solved reaction count. Output to feasibility.csv."

### Multi-objective route optimization
> "Use AiZynthFinder MCTS and only scorer names available in my installed configuration. Compare the top 10 routes by solved status, route length, and stock availability, reranking after search where necessary."

### Forward validation
> "For each route from AiZynthFinder, predict the forward reaction with Molecular Transformer; report whether the predicted product matches the target SMILES."

## What the Agent Will Do

1. Load AiZynthFinder and explicitly select the configured expansion-policy and stock keys used for the run.
2. Run MCTS with a recorded search configuration and budget.
3. Build routes from solved tree.
4. Compare routes by solved state, reaction count, in-stock leaves, and configured scores.
5. Optionally compare with a maintained or exactly version-pinned template-free/forward model.
6. Output route SMILES, building blocks, scorer values, reaction count, solved state, and stock provenance; report input-parse and execution failures separately from completed searches with no solved route.

## Tips

- Template coverage is domain-dependent; report the policy source/training domain and inspect failures rather than assuming a universal coverage percentage.
- Always check `in_stock` flag for each leaf; non-stock leaves require additional synthesis.
- AiZynthFinder uses MCTS; combine documented search rewards or rerank completed routes instead of selecting a nonexistent `mo_mcts` algorithm.
- For forward checks, report the observed top-k round-trip match rate for the selected model and reaction representation; do not assume a universal pass rate.
- For batch retrosynthesis on 1000+ compounds, use `aizynthcli` CLI.

## Related Skills

- chemoinformatics/molecular-io - Parse SMILES
- chemoinformatics/molecular-standardization - Standardize before retrosynthesis
- chemoinformatics/generative-design - Add feasibility to generative scoring
- chemoinformatics/reaction-enumeration - Forward direction
- chemoinformatics/admet-prediction - Filter targets before retrosynthesis

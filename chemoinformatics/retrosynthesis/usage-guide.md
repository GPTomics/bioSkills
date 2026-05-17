# Retrosynthesis Usage Guide

## Overview

Plan synthetic routes from target molecules back to commercial building blocks using AiZynthFinder 4.0 (template-based MCTS, AstraZeneca) or Chemformer (template-free Transformer). Score routes by depth, building-block availability, and forward-prediction validation.

## Prerequisites

```bash
pip install aizynthfinder chemformer
```

Download AiZynthFinder USPTO model + building-block stocks.

## Quick Start

Tell the AI agent what to do:
- "Plan retrosynthesis for this target SMILES using AiZynthFinder"
- "Batch retrosynthesize 100 generated molecules; report feasibility"
- "Find shorter routes using multi-objective MCTS"
- "Validate retrosynthesis with forward prediction via Molecular Transformer"
- "Estimate synthesis cost using building-block prices"

## Example Prompts

### Single target retrosynthesis
> "Plan retrosynthesis for SMILES 'CC(=O)Nc1ccc(C(=O)Nc2ccccc2)cc1' using AiZynthFinder. USPTO templates, ZINC stock. Output top 5 routes with depth, score, building-block list."

### Batch feasibility screening
> "Batch retrosynthesize 100 generated SMILES from candidates.csv. Report 'easy', 'feasible', or 'unsolved' per compound. Output to feasibility.csv."

### Multi-objective route optimization
> "Use MO-MCTS in AiZynthFinder to find routes balancing state_score (0.5), broken_bonds_score (0.3), and route_length (0.2, minimize). Report top 10 routes."

### Forward validation
> "For each route from AiZynthFinder, predict the forward reaction with Molecular Transformer; report whether the predicted product matches the target SMILES."

## What the Agent Will Do

1. Load AiZynthFinder with USPTO templates + ZINC building-block stock.
2. Run MCTS tree search (default 100 iterations, 120s time limit).
3. Build routes from solved tree.
4. Filter routes: depth, in-stock leaves, score.
5. Optionally validate with Chemformer or Molecular Transformer.
6. Output: route SMILES, building blocks, score, depth.

## Tips

- USPTO templates cover ~80% of common medchem; novel chemistry may need Chemformer.
- Always check `in_stock` flag for each leaf; non-stock leaves require additional synthesis.
- Multi-objective MCTS (MO-MCTS) avoids long suboptimal routes.
- Forward validation: ~30-50% of retrosynthesis routes pass round-trip check.
- For batch retrosynthesis on 1000+ compounds, use `aizynthcli` CLI.

## Related Skills

- chemoinformatics/molecular-io - Parse SMILES
- chemoinformatics/molecular-standardization - Standardize before retrosynthesis
- chemoinformatics/generative-design - Add feasibility to generative scoring
- chemoinformatics/reaction-enumeration - Forward direction
- chemoinformatics/admet-prediction - Filter targets before retrosynthesis

# Strain Design - Usage Guide

## Overview

Computational strain design searches a genome-scale model for gene or reaction interventions that make an organism overproduce a target chemical. The unifying idea is growth-coupling: rather than just deleting competing pathways (which the cell routes around and evolution erodes), a good design makes product secretion obligatory for growth, so selection in the bioreactor maintains production. OptKnock encodes this as a bilevel optimization - the inner problem is the cell maximizing its own growth, the outer is the engineer maximizing product at that optimum - but it is optimistic, assuming the cell cooperates; RobustKnock and minimal cut sets give more trustworthy, worst-case-aware designs. Every design is a hypothesis about a model with no regulation, kinetics, toxicity, or genetic stability, so it must be verified on the production envelope and then in the lab. The optimization is MILP and combinatorially hard, so bound the intervention set and use a strong solver.

## Prerequisites

```bash
pip install straindesign     # OptKnock/RobustKnock/MCS/OptCouple on COBRApy; MILP
pip install cameo            # heuristic (OptGene) and FSEOF over/under-expression targets (optional)
# MILP needs a solver: GLPK/SCIP (open source) work for small models; CPLEX/Gurobi (academic) for genome-scale.
```

Inputs: a curated genome-scale model, a target product exchange reaction, and a defined medium. Know the model's biomass reaction id (`linear_reaction_coefficients(model)`).

## Quick Start

Tell your AI agent:
- "Design knockouts to overproduce succinate in E. coli, coupled to growth"
- "Use RobustKnock instead of OptKnock so the design is robust to the cell's choices"
- "Compute minimal cut sets that force my product"
- "Find over-expression targets for my pathway with FSEOF"
- "Verify that my design is actually growth-coupled"

## Example Prompts

### Growth-Coupled Knockouts
> "On anaerobic glucose, find up to four knockouts that couple succinate secretion to growth in the E. coli core model, and verify with the production envelope that the strain cannot grow without making succinate."

### Robust vs Optimistic
> "My OptKnock design produces zero product at max growth. Explain OptKnock's optimism and redo it with RobustKnock so product is guaranteed in the worst-case internal state."

### Minimal Cut Sets
> "Enumerate the smallest sets of reaction knockouts that force my target product, and translate them to gene deletions via the GPR."

### Over-Expression
> "My bottleneck is low flux through an existing pathway, not a competing drain. Use FSEOF to find amplification targets rather than knockouts."

## What the Agent Will Do

1. Confirm the model is curated and identify the biomass and product reactions on a defined medium.
2. Build an SDModule (OptKnock/RobustKnock/MCS/OptCouple) with growth as the inner and product as the outer objective, plus a minimum-growth constraint.
3. Run `compute_strain_designs` with an intervention budget, solution cap, time limit, and a strong solver.
4. Decode reaction knockouts and translate them to gene deletions via the GPR, excluding essential genes.
5. Verify growth-coupling by pinning growth near its max and minimizing the product (nonzero minimum = coupled), or via the production envelope.
6. Prefer RobustKnock/MCS when the design must be trustworthy, and flag that the design is a model hypothesis to test in vivo.

## Tips

- Growth-coupling is the goal: a design where the cell cannot grow well without secreting product is what makes production evolutionarily stable.
- OptKnock is optimistic (it assumes the cell picks the growth-optimal state best for you); RobustKnock guards the worst case and MCS gives guarantees - prefer them when the design matters.
- Verify coupling directly: pin growth near max, minimize product; a nonzero minimum means production is obligatory. A design that lets product go to zero at max growth is weakly coupled.
- MILP is hard: set a `time_limit`, bound `max_cost` (interventions), cap `max_solutions`, and use CPLEX/Gurobi for genome-scale models - GLPK is slow.
- StrainDesign encodes a knockout as `reaction_id -> -1.0` in the design dict; translate reaction knockouts to gene deletions through the GPR before calling a design realizable.
- Do not design knockouts of essential genes; cross-check with gene-essentiality.
- Anaerobic/fermentative conditions make many products (succinate, ethanol, lactate) naturally easier to growth-couple than aerobic ones.
- FSEOF and cameo find over/under-expression targets - reach for these when the limit is low pathway flux rather than a competing drain a knockout would remove.
- A curated model is a prerequisite; a design on a model with an energy-generating cycle or wrong biomass is meaningless.
- The design is a hypothesis: FBA has no regulation, kinetics, toxicity, or genetic stability. Validate the envelope, then build and test the strain.

## Related Skills

- systems-biology/flux-balance-analysis - Production envelope and the FBA/medium foundation
- systems-biology/gene-essentiality - Avoid knocking out essential genes; GPR-to-gene mapping
- systems-biology/model-curation - A curated model is required for a trustworthy design
- systems-biology/context-specific-models - Constrain the chassis to a condition before designing
- metabolomics/pathway-mapping - Interpret the pathways a design affects

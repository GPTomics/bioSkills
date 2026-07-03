# Flux Balance Analysis - Usage Guide

## Overview

Flux balance analysis predicts metabolic phenotypes by solving a linear program over a genome-scale model: it maximizes a biomass (or custom) objective subject to steady-state mass balance and exchange bounds. The single most important thing to internalize is that FBA is under-determined - the maximal objective is reached across a whole face of the solution space, so `model.optimize()` returns one arbitrary flux vector among many. The objective value is robust; individual internal fluxes are not, and must be qualified with flux variability analysis (ranges), parsimonious FBA (a minimal-flux representative), or sampling (a distribution). A nonzero growth rate is only as trustworthy as the medium and the biomass reaction that produced it.

## Prerequisites

```bash
pip install cobra
# Optional faster/more robust solvers for large or MILP/QP work:
#   HiGHS ships with cobra; academic CPLEX or Gurobi are free and recommended for genome-scale models.
```

Inputs: a genome-scale model in SBML (`.xml`) or COBRA JSON, or a built-in test model (`load_model('textbook')`, `load_model('iJO1366')`). Curated models are at http://bigg.ucsd.edu/models.

## Quick Start

Tell your AI agent:
- "Run FBA on the E. coli core model and report the growth rate"
- "Set glucose minimal medium and predict growth, then tell me if it looks realistic"
- "Compute FVA ranges and flag which reactions are actually determined vs free"
- "Give me a parsimonious flux distribution for my model"
- "Sample the flux space and report the distribution for the key exchange reactions"

## Example Prompts

### Growth Prediction
> "Load iJO1366, set glucose-minimal aerobic medium with 10 mmol/gDW/h glucose, run FBA, and tell me whether the growth rate is biological or an artifact of an open exchange."

### Alternate Optima
> "Run FBA and then FVA at 100% and 90% of optimum on my model, and list the reactions whose flux is fully determined versus those free to vary - I want to know which single-FBA fluxes I can trust."

### Realistic Flux Distribution
> "Compare a plain FBA solution to pFBA and a loopless solution for my curated model, and explain which internal fluxes change and why."

### Flux Sampling
> "Uniformly sample the flux space of my model constrained to >=90% of max growth, and give me the median and interquartile range for the glycolysis and TCA reactions."

## What the Agent Will Do

1. Load the model (SBML/JSON/built-in) and confirm the objective is the intended biomass reaction.
2. Define the medium by closing all exchanges and opening only the intended uptakes (uptake = negative lower bound), or via the `model.medium` dict.
3. Solve the LP with `optimize()` / `slim_optimize()` and report the objective value with its solver status.
4. Qualify internal fluxes with FVA (ranges), pFBA (parsimonious representative), or sampling (distributions) rather than trusting a single flux vector.
5. Apply loopless correction where internal cycles inflate flux magnitudes.
6. Flag artifacts: open-exchange-inflated growth, blocked reactions, solver-tolerance noise.

## Tips

- Always state the medium and biomass reaction alongside any growth number; the value is meaningless without them.
- Use `slim_optimize()` inside loops (medium scans, comparisons) - it returns only the objective float and is much faster than building a full `Solution`.
- Wrap temporary bound changes in `with model:` so they revert automatically and comparisons never leak state.
- Uptake is a negative flux by convention; open an uptake by setting the exchange's lower bound negative, not its upper bound.
- Do not report `solution.fluxes[rxn]` for an internal reaction without checking its FVA range - a wide range means the value was arbitrary.
- pFBA gives a single tidy representative but does not make the fluxes "true"; sampling is the honest choice when uncertainty matters.
- Check sampler mixing (multiple seeds/chains should agree) before trusting a sampled distribution; short chains are autocorrelated and misleading.
- Plain FBA cannot reproduce overflow metabolism (acetate at high glucose, Crabtree) - that needs an enzyme/proteome-constrained model (GECKO, sMOMENT).
- FBA predicts yields and growth well but internal flux magnitudes poorly; validate real intracellular fluxes with 13C-MFA (metabolomics/isotope-tracing), not FBA alone.
- Prefer HiGHS or academic CPLEX/Gurobi over GLPK for genome-scale models and any loopless/MOMA/MILP work; GLPK is the most likely to return degenerate or marginally-infeasible answers.

## Related Skills

- systems-biology/gene-essentiality - Knockout screens; MOMA/ROOM for immediate (non-re-optimized) mutants
- systems-biology/context-specific-models - Constrain FBA with tissue/condition expression data
- systems-biology/strain-design - Turn the production envelope into knockout designs
- systems-biology/community-metabolic-modeling - FBA over multi-species communities
- metabolomics/isotope-tracing - Measure intracellular fluxes with 13C-MFA

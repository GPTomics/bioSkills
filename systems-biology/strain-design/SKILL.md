---
name: bio-systems-biology-strain-design
description: Computes metabolic-engineering strain designs on genome-scale models with StrainDesign (OptKnock, RobustKnock, minimal cut sets, OptCouple) and cameo (heuristic knockout and FSEOF over/under-expression targets), finding gene/reaction interventions that couple product formation to growth. Use when designing knockouts to overproduce a target chemical, choosing between OptKnock and RobustKnock, growth-coupling a product so evolution maintains it, computing minimal cut sets, finding amplification targets with FSEOF, or understanding why MILP strain design needs a strong solver and why a design is only a hypothesis.
tool_type: python
primary_tool: straindesign
---

## Version Compatibility

Reference examples tested with: StrainDesign 1.15+, COBRApy 0.29+, Python 3.10+ (cameo optional)

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Note: OptKnock/RobustKnock/MCS are MILP problems and are far harder than plain FBA; StrainDesign supports GLPK/SCIP (open source) and CPLEX/Gurobi (academic, much faster and more robust for genome-scale). Set a `time_limit`. Reaction-based designs must be translated back to gene knockouts via GPRs.

# Strain Design

**"Design knockouts to make my organism overproduce a chemical"** -> Search for a set of gene/reaction interventions that couples product formation to growth, so the engineered strain cannot grow well without secreting the target.
- Python: `straindesign.compute_strain_designs(model, sd_modules=[SDModule(model, OPTKNOCK, ...)])`; cameo for heuristics/FSEOF

## The governing principle: growth-coupling makes evolution enforce the design

The core idea of computational strain design is growth-coupling. A naive "just delete the competing pathways" design is fragile: the cell will find an alternate flux route, or evolution in the bioreactor will erode production because making product costs the cell resources. A growth-COUPLED design instead makes product secretion obligatory for growth - the cell physically cannot reach high growth without also secreting the target, so selection maintains production instead of eroding it. This is why OptKnock is a BILEVEL optimization: the inner problem is the cell maximizing its own growth, the outer problem is the engineer maximizing product AT that inner optimum. Consequences:

- OptKnock is optimistic: it assumes the cell, among its growth-optimal states, picks the one best for the engineer. The cell need not. RobustKnock fixes this by maximizing product in the WORST-case inner optimum - a more conservative, more trustworthy design.
- A design is a HYPOTHESIS about a model, not a strain. FBA has no regulation, no enzyme kinetics, no toxicity, no genetic stability; a computationally growth-coupled design can fail in construction or in the bioreactor. Validate with the production envelope, then in the lab.
- MILP strain design is combinatorially hard. Bound the intervention set (`max_cost`), cap solutions, set a time limit, and use a strong solver (CPLEX/Gurobi for genome-scale). Reaction knockouts must be mapped back to gene deletions through the GPR to be realizable.

## Decision: which strain-design method

| Goal | Method | Trade-off |
|------|--------|-----------|
| Growth-coupled knockouts, optimistic | OptKnock (Burgard 2003) | bilevel; assumes the cell cooperates at its growth optimum |
| Growth-coupled knockouts, conservative | RobustKnock (Tepper & Shlomi 2010) | guarantees product in the worst-case inner optimum; harder |
| Guaranteed intervention sets, enumerate all minimal | Minimal Cut Sets (von Kamp & Klamt 2014) | strong guarantees; enumerates smallest intervention sets |
| Strong growth-coupling (obligatory) | OptCouple | maximizes the growth-coupling potential directly |
| Over/under-EXPRESSION targets, not just knockouts | FSEOF (Choi 2010) / cameo | scans fluxes that rise with enforced product; amplification targets |
| Heuristic/evolutionary search when MILP is intractable | OptGene / cameo | fast approximate designs; no optimality guarantee |

Prefer RobustKnock or MCS over plain OptKnock when the design must be trustworthy; OptKnock's optimism is a well-known way to overstate a design.

## Growth-Coupled Knockouts with StrainDesign (OptKnock)

**Goal:** Find a small set of reaction knockouts that couples secretion of a target product to growth.

**Approach:** Build an OptKnock `SDModule` with the cell's growth as the inner objective and product secretion as the outer objective, plus a minimum-growth constraint so the design keeps the strain viable, then call `compute_strain_designs` with an intervention budget and solver. Translate the returned reaction knockouts back to gene deletions via the GPR.

```python
import cobra
import straindesign as sd

model = cobra.io.load_model('textbook')
biomass = 'Biomass_Ecoli_core'   # the model's actual biomass reaction id (verify per model)

optknock = sd.SDModule(
    model, sd.OPTKNOCK,
    inner_objective=biomass,            # the cell maximizes growth
    outer_objective='EX_ac_e',          # the engineer maximizes acetate secretion
    constraints=[f'{biomass} >= 0.3'],  # keep the strain viable
)

solutions = sd.compute_strain_designs(
    model, sd_modules=[optknock],
    max_cost=3,          # at most 3 interventions
    max_solutions=3,
    solver='glpk',       # use 'cplex'/'gurobi' for genome-scale models
    time_limit=120,
)
# solutions.reaction_sd is a list of intervention dicts {reaction_id: marker}; a knockout is
# marked -1.0 (not 0). Verify this marker for the installed StrainDesign version -- a wrong marker
# silently yields empty designs. For a knockout-only OptKnock module every entry is a knockout.
for design in solutions.reaction_sd:
    print('knockouts:', [rid for rid, mark in design.items() if mark == -1.0])
```

## Verify Growth-Coupling with the Production Envelope

```python
from cobra.flux_analysis import production_envelope

# A genuinely growth-coupled design shows a NONZERO minimum product flux across the growth range:
# the strain cannot grow without secreting product. Apply the design's knockouts, then:
env = production_envelope(model, reactions=['EX_ac_e'])   # objective defaults to biomass
# Inspect the lower bound of product at each growth level; if it can be zero at max growth, the
# coupling is weak (the OptKnock-optimism problem) -- consider RobustKnock. See flux-balance-analysis.
```

## Over-Expression Targets (FSEOF, cameo)

```python
# Knockouts are not the only lever. FSEOF (flux scanning with enforced objective flux) finds
# reactions whose flux RISES as product formation is enforced -- candidate amplification/over-
# expression targets. cameo implements FSEOF and heuristic (evolutionary) design search:
#   from cameo.strain_design import OptGene           # heuristic knockout search
#   from cameo.strain_design.deterministic import FSEOF
# Use FSEOF/over-expression when the bottleneck is low flux through an existing pathway rather than
# a competing drain that a knockout would remove.
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| `compute_strain_designs` never finishes | MILP is hard and GLPK is slow on genome-scale | set `time_limit`, lower `max_cost`, use CPLEX/Gurobi |
| Design gives zero product when built | OptKnock optimism: the cell chose a different growth-optimal state | use RobustKnock, or check the production envelope's lower bound |
| Constraint parser rejects the biomass id | wrong reaction id string for this model | look up the actual objective reaction id (`linear_reaction_coefficients`) |
| Design not realizable in the lab | reaction knockouts have no clean gene mapping, or hit an essential gene | translate reaction KOs to gene KOs via GPR; exclude essential genes |
| Predicted overproduction never materializes | FBA has no regulation/kinetics/toxicity/stability | treat the design as a hypothesis; validate the envelope, then in vivo |
| No feasible design found | growth constraint too tight or product infeasible on the medium | relax the minimum-growth constraint; confirm the product can be made on the medium |

## Related Skills

- systems-biology/flux-balance-analysis - Production envelope and the FBA/medium foundation
- systems-biology/gene-essentiality - Avoid designing knockouts of essential genes; GPR mapping
- systems-biology/model-curation - A curated model is a prerequisite for a trustworthy design
- systems-biology/context-specific-models - Constrain the chassis to a condition before designing
- metabolomics/pathway-mapping - Interpret the affected pathways of a design

## References

- Burgard AP, Pharkya P, Maranas CD. 2003. OptKnock: a bilevel programming framework for identifying gene knockout strategies for microbial strain optimization. *Biotechnol Bioeng* 84(6):647-657.
- Schneider P, Bekiaris PS, von Kamp A, Klamt S. 2022. StrainDesign: a comprehensive Python package for computational design of metabolic networks. *Bioinformatics* 38(21):4981-4983.
- Tepper N, Shlomi T. 2010. Predicting metabolic engineering knockout strategies for chemical production: accounting for competing pathways. *Bioinformatics* 26(4):536-543. (RobustKnock)
- von Kamp A, Klamt S. 2014. Enumeration of smallest intervention strategies in genome-scale metabolic networks. *PLoS Comput Biol* 10(1):e1003378. (minimal cut sets)
- Choi HS, Lee SY, Kim TY, Woo HM. 2010. In silico identification of gene amplification targets for improvement of lycopene production. *Appl Environ Microbiol* 76(10):3097-3105. (FSEOF)
- Cardoso JGR, Jensen K, Lieven C, et al. 2018. Cameo: a Python library for computer-aided metabolic engineering and optimization of cell factories. *ACS Synth Biol* 7(4):1163-1166.

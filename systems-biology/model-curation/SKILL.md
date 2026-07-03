---
name: bio-systems-biology-model-curation
description: Validates, gap-fills, and standardizes genome-scale metabolic models using memote for consistency and annotation scoring and COBRApy for manual curation, including mass/charge balance, energy-generating-cycle detection, dead-end resolution, GPR fixes, and SBML/SBO/MIRIAM annotation. Use when improving a draft model, gap-filling to a target medium, detecting erroneous ATP-from-nothing cycles, interpreting a memote score correctly (consistency vs biological validity), validating predictions against measured growth/essentiality, or preparing a model for publication.
tool_type: python
primary_tool: memote
---

## Version Compatibility

Reference examples tested with: memote 0.17+, COBRApy 0.29+, Python 3.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Note: a memote SCORE is comparable only within a memote version (the test suite evolves). The Python entry points are `memote.suite.api.test_model` / `snapshot_report` (not `run`/`snapshot`), and `cobra.flux_analysis.gapfill` has no `demand` argument (it gap-fills toward the model's objective).

# Model Curation

**"Validate and improve the quality of my metabolic model"** -> Score the model for consistency and annotation with memote, then use COBRApy to fix mass/charge imbalances, remove energy-generating cycles, gap-fill deliberately, and validate predictions against data.
- CLI: `memote report snapshot model.xml --filename report.html`
- Python: `cobra.flux_analysis.gapfill()`, mass/charge balance and energy-cycle checks (COBRApy)

## The governing principle: memote scores consistency, not biological correctness

The most common misconception in the field is that a high memote score means a good model. memote's scored total is a weighted sum of stoichiometric consistency, mass/charge balance, annotation coverage (KEGG/ChEBI/BiGG), SBO-term presence, and SBML/FBC conformance. All of that is HYGIENE - it measures whether the model is well-formed and well-annotated, NOT whether it predicts biology. A model can score 90% and mispredict every knockout and every growth phenotype. Optimizing the score for its own sake is Goodhart's law made concrete.

Curation therefore has two separable axes, and both must be reported:

1. Syntactic/consistency (what memote measures): mass and charge balance, no stoichiometric leaks, annotation, SBO terms, no blocked reactions, and - the one predictive-adjacent test - no erroneous energy-generating cycles.
2. Predictive validity (what memote does NOT measure): does the model reproduce measured growth rates, carbon-source usage (Biolog), and gene essentiality on the matched medium? This is the step that separates a curated model from a merely tidy one, and it lives outside memote.

The single most dangerous defect a high score can hide is an energy-generating cycle: a set of reactions that produces ATP/NADH from nothing. It arises from reversible reactions and blind gap-filling, passes mass balance, inflates growth, and invalidates every flux prediction. Test for it explicitly.

## Decision: which curation action for which symptom

| Symptom | Action | Tool |
|---------|--------|------|
| Model cannot grow on the target medium | gap-fill toward the objective, from a universal DB | `cobra.flux_analysis.gapfill` |
| Growth is implausibly high / ATP from nothing | detect and break energy-generating cycles | max-ATP-with-no-uptake test; constrain directionality |
| Reactions unbalanced | fix mass/charge (usually protons at pH 7) | per-reaction element/charge sum |
| Metabolite never produced or never consumed | resolve dead-end (add reaction or fix stoichiometry) | connectivity scan |
| Low annotation / SBO score | add MIRIAM annotations and SBO terms | memote report + manual/annotation tools |
| Predictions wrong despite high score | validate against measured growth/essentiality | separate experimental comparison (not memote) |

## memote: score consistency, then read the report, not just the number

```bash
pip install memote

memote run model.xml                                              # run the test suite (pytest-based)
memote report snapshot model.xml --filename report.html          # human-readable HTML report
```

```python
# Programmatic entry points (verify against the installed memote version):
from memote.suite.api import test_model, snapshot_report

code, result = test_model(model, results=True)   # result is a MemoteResult (the raw test outcomes)
html = snapshot_report(result, html=True)         # render the same report programmatically
# Read WHICH tests fail (consistency, energy cycles, unbalanced reactions) -- the total % is not
# a measure of biological correctness.
```

## Detect Energy-Generating Cycles (the defect a high score hides)

**Goal:** Prove the model cannot manufacture any energy currency (ATP, NADH, NADPH, FADH2, ...) from nothing.

**Approach:** Close every exchange so no nutrients enter and zero the ATP-maintenance lower bound (its NGAM floor would otherwise make a closed model infeasible for the wrong reason). Then, for EACH energy currency, add a moiety-conserving dissipation reaction (charged -> discharged) and maximize it. A result of 0 (or infeasible) per currency is correct; any positive finite flux is an erroneous energy-generating cycle to trace and fix by constraining reaction directionality. EGCs are not ATP-only, so the sweep must cover every currency present (Fritzemeier 2017); proton-motive-force cycles are subtler and are handled by memote's dedicated EGC test.

```python
# Dissipation stoichiometry per currency (BiGG ids); genome-scale models also test GTP/CTP/UTP/q8h2.
DISSIPATION = {'atp': {'atp_c': -1, 'h2o_c': -1, 'adp_c': 1, 'pi_c': 1, 'h_c': 1},
               'nadh': {'nadh_c': -1, 'nad_c': 1, 'h_c': 1},
               'nadph': {'nadph_c': -1, 'nadp_c': 1, 'h_c': 1}}

def energy_generating_cycles(model, dissipations=DISSIPATION):
    '''Max free-charging flux per currency with ALL uptake closed; >0 => energy-generating cycle.'''
    out = {}
    with model:
        for ex in model.exchanges:
            ex.lower_bound = 0                                   # no nutrients at all
        if 'ATPM' in model.reactions:
            model.reactions.get_by_id('ATPM').lower_bound = 0    # remove the NGAM floor before testing
        for name, stoich in dissipations.items():
            if any(m not in model.metabolites for m in stoich):
                continue                                         # currency absent from this model
            with model:
                r = cobra.Reaction(f'EGC_{name}')
                r.add_metabolites({model.metabolites.get_by_id(m): c for m, c in stoich.items()})
                r.bounds = (0, 1000)
                model.add_reactions([r])
                model.objective = r
                out[name] = model.slim_optimize()   # 0/infeasible = OK; positive finite = cycle
    return out
```

## Gap-Fill Toward the Objective (deliberately, on a stated medium)

**Goal:** Add the fewest reactions from a universal database that let the model grow on a defined medium.

**Approach:** Set the medium and the biomass objective, then call `gapfill` (which minimizes added reactions to reach the objective at `lower_bound`). There is no `demand` argument; `demand_reactions=False` avoids adding demand reactions for every metabolite. Record and low-confidence-flag every added reaction.

```python
from cobra.flux_analysis import gapfill

universal = cobra.io.read_sbml_model('universal_model.xml')   # e.g. a BiGG universal model
solutions = gapfill(model, universal, lower_bound=0.05, demand_reactions=False, iterations=3)
for i, rxns in enumerate(solutions):
    print(f'solution {i+1}: {[r.id for r in rxns]}')          # alternative gap-fill sets
# Adding these forces growth; that is not evidence they are biologically present. Flag them.
```

## Mass and Charge Balance

```python
def imbalance(reaction):
    '''Return the element and charge imbalance of a reaction (empty dict + 0 charge if balanced).'''
    mass = {}
    charge = 0
    for met, coef in reaction.metabolites.items():
        if met.formula:
            for element, n in met.elements.items():
                mass[element] = mass.get(element, 0) + coef * n
        if met.charge is not None:
            charge += coef * met.charge
    return {e: v for e, v in mass.items() if abs(v) > 1e-6}, charge

# This is a PER-REACTION element/charge check. It is distinct from stoichiometric CONSISTENCY --
# a whole-network LP (Gevorgyan 2008, what memote tests) that finds mass leaks without needing
# formulas. "All reactions mass-balanced" does not imply the network is stoichiometrically consistent.
# Exchange/demand/sink AND the biomass pseudo-reaction are intentionally imbalanced; skip them.
# Proton (H) imbalance at pH 7 is the most common real fix.
from cobra.util.solver import linear_reaction_coefficients
pseudo = set(model.boundary) | set(linear_reaction_coefficients(model))   # boundary + objective (biomass)
unbalanced = [(r.id, imbalance(r)) for r in model.reactions
              if r not in pseudo and (imbalance(r)[0] or abs(imbalance(r)[1]) > 1e-6)]
```

## Dead-End Metabolites

```python
def dead_ends(model):
    '''Metabolites that can only be produced or only consumed (a network gap or wrong stoichiometry).'''
    out = []
    for met in model.metabolites:
        produced = any(r.metabolites[met] > 0 for r in met.reactions)
        consumed = any(r.metabolites[met] < 0 for r in met.reactions)
        if not (produced and consumed):
            out.append(met.id)
    return out
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| "High memote score, so the model is good" | score measures consistency/annotation, not prediction | validate against measured growth/essentiality separately |
| Growth is huge; ATP looks free | erroneous energy-generating cycle | run the max-ATP-with-no-uptake test; constrain the offending reactions' directionality |
| `gapfill(... demand=...)` TypeError | there is no `demand` argument | set the objective and use `lower_bound=`/`demand_reactions=False` |
| `memote.suite.api.run`/`snapshot` AttributeError | wrong names | use `test_model` / `snapshot_report` |
| Many reactions flagged unbalanced | protons/charge at pH 7, or exchange reactions counted | skip exchange/sink/demand; fix H and charge first |
| Model still mispredicts after high score | consistency fixed, biology not validated | compare to Biolog carbon sources and an essentiality screen on the matched medium |

## Related Skills

- systems-biology/metabolic-reconstruction - Produces the draft this skill curates
- systems-biology/flux-balance-analysis - Test the curated model's predictions
- systems-biology/gene-essentiality - Validate curation against measured essentiality
- pathway-analysis/kegg-pathways - Source KEGG annotations for reactions/metabolites
- database-access/uniprot-access - Cross-reference gene/protein annotations

## References

- Lieven C, Beber ME, Olivier BG, et al. 2020. MEMOTE for standardized genome-scale metabolic model testing. *Nat Biotechnol* 38(3):272-276.
- Fritzemeier CJ, Hartleb D, Szappanos B, Papp B, Lercher MJ. 2017. Erroneous energy-generating cycles in published genome-scale metabolic networks: identification and removal. *PLoS Comput Biol* 13(4):e1005494.
- Noor E, Haraldsdottir HS, Milo R, Fleming RMT. 2013. Consistent estimation of Gibbs energy using component contributions. *PLoS Comput Biol* 9(7):e1003098. (thermodynamic directionality)
- Thiele I, Palsson BO. 2010. A protocol for generating a high-quality genome-scale metabolic reconstruction. *Nat Protoc* 5(1):93-121.
- Orth JD, Palsson BO. 2010. Systematizing the generation of missing metabolic knowledge. *Biotechnol Bioeng* 107(3):403-412. (gap analysis)
- Ebrahim A, Lerman JA, Palsson BO, Hyduke DR. 2013. COBRApy: constraint-based reconstruction and analysis for Python. *BMC Syst Biol* 7:74.

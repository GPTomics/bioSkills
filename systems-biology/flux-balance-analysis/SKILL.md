---
name: bio-systems-biology-flux-balance-analysis
description: Performs flux balance analysis (FBA), flux variability analysis (FVA), parsimonious FBA (pFBA), loopless FBA, flux sampling, and production envelopes on genome-scale metabolic models with COBRApy, solving the biomass-maximization linear program under a defined medium. Use when predicting growth rate on a carbon source, computing flux ranges and alternative optima (FVA), setting exchange bounds and minimal media, distinguishing a real growth phenotype from an under-constrained model, sampling the flux solution space, or choosing between FBA, pFBA, loopless FBA, and sampling for a flux distribution.
tool_type: python
primary_tool: cobrapy
---

## Version Compatibility

Reference examples tested with: COBRApy 0.29+, Python 3.10+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Note: the objective LP is solved by a backend solver. GLPK (the COBRApy default) can return degenerate or marginally-infeasible answers on large or ill-conditioned models; prefer HiGHS (if available in the installed cobra/optlang build) or a free academic CPLEX/Gurobi for genome-scale work and any MILP/QP method (loopless, MOMA). Set with `model.solver = 'glpk'|'highs'|'gurobi'|'cplex'`.

# Flux Balance Analysis

**"Predict growth rate and metabolic fluxes for my organism"** -> Solve a linear program over a genome-scale metabolic model that maximizes a biomass (or custom) objective subject to steady-state mass balance and exchange bounds, then quantify how much of that solution is actually determined.
- Python: `model.optimize()`, `cobra.flux_analysis.flux_variability_analysis()`, `pfba()`, `cobra.sampling.sample()` (COBRApy)

## The governing principle: FBA is an under-determined LP, so one solution is not "the" flux distribution

FBA imposes steady state (S*v = 0) plus bounds and maximizes an objective. The steady-state constraint is a pseudo-steady-state on fast-turnover metabolite POOLS, not on the organism. Critically, the optimal objective is usually reached on a whole FACE of the solution polytope, not a single point: many different internal flux vectors give the identical maximal growth. `model.optimize()` returns ONE arbitrary vertex of that face. Consequences that govern every downstream decision:

- The biomass/objective VALUE is typically robust and reproducible; individual internal FLUXES are often NOT. Never report a single `solution.fluxes[rxn]` as "the" flux without FVA (the range) or pFBA (a parsimonious representative) or sampling (the distribution).
- A nonzero growth rate is only as meaningful as the MEDIUM and the biomass reaction. An open or under-constrained exchange set inflates growth; a copied/generic biomass reaction encodes another organism's composition. "It grows" is often an artifact of an open exchange, not biology (see Common Errors).
- FBA predicts YIELDS and growth/essentiality well but internal flux magnitudes poorly. To validate actual intracellular fluxes, 13C metabolic flux analysis MEASURES them (metabolomics/isotope-tracing) - FBA does not.

## Decision: which analysis for which question

| Question | Method | Why |
|----------|--------|-----|
| Max growth rate / yield on a medium | FBA `model.optimize()` | single LP; objective value is the robust output |
| Is a reaction's flux determined, or free to vary? | FVA `flux_variability_analysis` | reports min/max flux at (near-)optimal growth; exposes alternate optima |
| One realistic representative flux vector | pFBA `pfba` | among optima, the one minimizing total flux (proxy for minimal enzyme cost) |
| Remove thermodynamically infeasible internal cycles | loopless (`loopless_solution`, or FVA `loopless=True`) | strips net flux around closed loops with no driving force |
| Full uncertainty / flux DISTRIBUTIONS, no objective needed | sampling `cobra.sampling.sample` | uniformly samples the solution space; use when the objective is unknown or confidence intervals are needed |
| Growth-vs-byproduct tradeoff for engineering | `production_envelope` | Pareto frontier of growth vs product secretion |
| Immediate knockout mutant flux (not re-optimized) | MOMA/ROOM -> gene-essentiality | minimal-adjustment, not re-optimization; see systems-biology/gene-essentiality |
| Overflow metabolism (acetate/Crabtree) missing | enzyme-constrained model (GECKO/sMOMENT) | plain FBA has no proteome budget; needs capacity constraints |

## Load Models

```python
import cobra

model = cobra.io.load_model('textbook')   # E. coli core (e_coli_core), 95 reactions, ships with cobra
model = cobra.io.load_model('iJO1366')    # genome-scale E. coli, 2583 reactions

model = cobra.io.read_sbml_model('model.xml')   # SBML (the standard exchange format)
model = cobra.io.load_json_model('model.json')  # COBRA JSON

# Curated genome-scale models: http://bigg.ucsd.edu/models (King 2016). Record the model
# version; predictions are only comparable within the same model release.
```

## Basic FBA and honest growth interpretation

**Goal:** Predict the maximum growth rate and a flux distribution under a defined medium, and judge whether the number is biological.

**Approach:** Load a model, confirm the objective is the intended biomass reaction, solve the LP, then read the objective value while treating individual fluxes as provisional until FVA/sampling confirms them.

```python
import cobra

model = cobra.io.load_model('textbook')

solution = model.optimize()
print(f'Objective (growth): {solution.objective_value:.4f} /h  status: {solution.status}')
print('Objective reaction:', str(model.objective.expression).split('*')[1].split()[0])

# A growth rate is interpretable ONLY against a stated medium and biomass reaction.
# Compare to a measured doubling time (mu = ln2 / t_double) rather than to a fixed
# "fast/slow" scale; absolute values are organism- and biomass-definition-specific.
```

## Set Medium (exchange bounds; uptake is a NEGATIVE lower bound)

**Goal:** Impose a defined nutrient environment so growth reflects the intended condition, not leftover open exchanges.

**Approach:** Close every exchange, then open only the intended uptakes. By COBRApy convention an exchange `EX_x_e` has flux < 0 for uptake and > 0 for secretion, so uptake is set through the lower bound. The `model.medium` dict is the concise idiom (its values are uptake magnitudes, positive).

```python
def set_minimal_medium(model, carbon_source='EX_glc__D_e', carbon_uptake=10):
    '''Close all uptake, then open a defined minimal medium.

    carbon_uptake in mmol/gDW/h; 10 is the standard E. coli aerobic glucose rate (iJO1366).
    '''
    for rxn in model.exchanges:
        rxn.lower_bound = 0

    minimal = {'EX_o2_e': 1000, 'EX_h2o_e': 1000, 'EX_h_e': 1000, 'EX_nh4_e': 1000,
               'EX_pi_e': 1000, 'EX_so4_e': 1000, 'EX_k_e': 1000, 'EX_mg2_e': 1000}
    for ex_id, uptake in minimal.items():
        if ex_id in model.reactions:
            model.reactions.get_by_id(ex_id).lower_bound = -uptake
    if carbon_source in model.reactions:
        model.reactions.get_by_id(carbon_source).lower_bound = -carbon_uptake
    return model

# The with-block reverts all bound changes on exit, so comparisons never leak state.
for cs in ['EX_glc__D_e', 'EX_ac_e', 'EX_succ_e']:
    with model:
        set_minimal_medium(model, carbon_source=cs)
        print(f'{cs}: growth = {model.slim_optimize():.4f}')   # slim_optimize returns the objective float only (fast)
```

## Flux Variability Analysis (FVA): expose alternate optima

```python
from cobra.flux_analysis import flux_variability_analysis

# Range each reaction can carry while holding the objective at (fraction_of_optimum) of the max.
fva = flux_variability_analysis(model, fraction_of_optimum=1.0)

# fraction_of_optimum < 1 relaxes the objective and reveals the alternative-optima span.
fva90 = flux_variability_analysis(model, fraction_of_optimum=0.9)

# loopless=True removes thermodynamically infeasible internal loops from the ranges (slower).
fva_ll = flux_variability_analysis(model, loopless=True)

# A reaction with maximum == minimum is fully determined; a wide range means the single
# FBA value for it was arbitrary. Blocked reactions have min == max == 0.
fva['range'] = fva['maximum'] - fva['minimum']
fva['blocked'] = fva['range'].abs() < 1e-9
```

## Parsimonious FBA (pFBA): one realistic representative

```python
from cobra.flux_analysis import pfba

# Among all optima, pFBA returns the flux vector minimizing the sum of absolute fluxes,
# an Occam's-razor proxy for minimal total enzyme cost (Lewis 2010). It is a principled
# single representative of the optimal face; it is NOT more "true" than the face itself,
# and it hugs the polytope boundary (a min-flux vertex), so report the FVA range alongside it
# when the internal flux values matter.
pfba_solution = pfba(model)
print(f'FBA total flux : {model.optimize().fluxes.abs().sum():.1f}')
print(f'pFBA total flux: {pfba_solution.fluxes.abs().sum():.1f}')
```

## Loopless FBA

```python
from cobra.flux_analysis import loopless_solution

# Projects an FBA solution onto a loopless one: no net flux around a closed cycle that
# lacks a thermodynamic driving force (Schellenberger 2011). Internal cycles are a common
# artifact of reversible reactions and gap-filling and inflate apparent flux magnitudes.
loopless = loopless_solution(model)
```

## Flux Sampling: the distribution, without choosing an objective

**Goal:** Characterize the whole space of feasible steady-state fluxes rather than one optimum, giving each reaction a distribution and confidence interval.

**Approach:** Uniformly sample the (optionally objective-constrained) solution polytope with a Markov-chain sampler (OptGP or ACHR). Use when there is no clear objective, when the objective face is large, or when uncertainty on internal fluxes matters more than an optimum.

```python
from cobra.sampling import sample

# n samples; method 'optgp' (parallel) or 'achr'. Thinning reduces autocorrelation.
samples = sample(model, n=1000, method='optgp', thinning=100, seed=1)
print(samples['PFK'].describe())   # per-reaction distribution, e.g. median and IQR

# To sample only high-growth states, constrain the objective first (e.g. biomass >= 0.9*max)
# inside a `with model:` block, then sample. Check mixing before trusting the distribution
# (multiple chains/seeds should agree); short chains give correlated, misleading samples.
```

## Production Envelope (growth vs product tradeoff)

```python
from cobra.flux_analysis import production_envelope

# Pareto frontier of growth against a secreted product; the design space for strain engineering.
# objective defaults to the model's objective (the biomass reaction) when omitted.
env = production_envelope(model, reactions=['EX_ac_e'])
# To actually DESIGN knockouts that couple product to growth, see systems-biology/strain-design.
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Growth is 0 / status 'infeasible' | medium closed, or an essential exchange left at lb=0, or a biomass precursor unproducible | check `model.medium`; open the minimal exchange set; confirm biomass precursors have a route |
| Suspiciously high growth on "minimal" media | an exchange left open to uptake (rich carbon, or all exchanges default-open) | close all exchange lower bounds to 0, then open only the intended set; audit `model.medium` |
| A reaction's flux changes every run / disagrees between tools | alternate optima - the value was one arbitrary vertex | report the FVA range or a pFBA/sampling value, not a single `optimize()` flux |
| Implausibly large internal fluxes | thermodynamically infeasible internal cycle | use `loopless_solution` or FVA `loopless=True` |
| Model predicts no acetate overflow at high glucose | plain FBA has no proteome/enzyme budget | use an enzyme-constrained model (GECKO/sMOMENT); FBA cannot see overflow |
| `phenotype_phase_plane(...)` raises TypeError | it is now a module, not a callable, in modern COBRApy | use `production_envelope` instead |
| Solver returns tiny nonzero "fluxes" that should be 0 | GLPK feasibility tolerance / degeneracy | switch to HiGHS or CPLEX/Gurobi; threshold fluxes at ~1e-6 |

## Related Skills

- systems-biology/gene-essentiality - Knockout screens, MOMA/ROOM for non-re-optimized mutants
- systems-biology/context-specific-models - Constrain FBA with expression data
- systems-biology/strain-design - Design knockouts from the production envelope
- systems-biology/community-metabolic-modeling - FBA over multi-species communities
- metabolomics/isotope-tracing - 13C-MFA to measure (not predict) intracellular fluxes
- metabolomics/pathway-mapping - Map measured metabolites onto model reactions

## References

- Orth JD, Thiele I, Palsson BO. 2010. What is flux balance analysis? *Nat Biotechnol* 28(3):245-248.
- Ebrahim A, Lerman JA, Palsson BO, Hyduke DR. 2013. COBRApy: constraint-based reconstruction and analysis for Python. *BMC Syst Biol* 7:74.
- Mahadevan R, Schilling CH. 2003. The effects of alternate optimal solutions in constraint-based genome-scale metabolic models. *Metab Eng* 5(4):264-276.
- Lewis NE, Hixson KK, Conrad TM, et al. 2010. Omic data from evolved E. coli are consistent with computed optimal growth from genome-scale models. *Mol Syst Biol* 6:390. (pFBA)
- Schellenberger J, Lewis NE, Palsson BO. 2011. Elimination of thermodynamically infeasible loops in steady-state metabolic models. *Biophys J* 100(3):544-553. (loopless)
- Megchelenbrink W, Huynen M, Marchiori E. 2014. optGpSampler: an improved tool for uniformly sampling the solution-space of genome-scale metabolic networks. *PLoS One* 9(2):e86587.
- Schuetz R, Kuepfer L, Sauer U. 2007. Systematic evaluation of objective functions for predicting intracellular fluxes in E. coli. *Mol Syst Biol* 3:119. (objective choice)
- Ibarra RU, Edwards JS, Palsson BO. 2002. Escherichia coli K-12 undergoes adaptive evolution to achieve in silico predicted optimal growth. *Nature* 420(6912):186-189.
- King ZA, Lu JS, Drager A, et al. 2016. BiGG Models: a platform for integrating, standardizing and sharing genome-scale models. *Nucleic Acids Res* 44(D1):D515-D522.

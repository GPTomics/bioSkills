---
name: bio-systems-biology-community-metabolic-modeling
description: Builds and simulates multi-species metabolic community models from member genome-scale models, using MICOM for abundance-weighted steady-state community FBA and cooperative tradeoff, SMETANA for cross-feeding and competition scoring, and SteadyCom/COMETS for common-growth-rate and dynamic simulation. Use when modeling a microbiome or co-culture, predicting cross-feeding and competition, abundance-weighting members from metagenomics, choosing steady-state vs dynamic community modeling, avoiding the compartment-pooling artifact, or judging how member-model quality and namespace propagate into community predictions.
tool_type: python
primary_tool: micom
---

## Version Compatibility

Reference examples tested with: MICOM 0.33+, COBRApy 0.29+, Python 3.10+ (SMETANA and COMETS are separate installs)

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Note: community FBA needs a QP solver for MICOM's cooperative tradeoff (HiGHS/CPLEX/Gurobi). Member models must share a namespace (BiGG vs ModelSEED), reconciled via MetaNetX before combining. SMETANA is a separate CLI (github.com/cdanielmachado/smetana); COMETS uses the cometspy toolbox.

# Community Metabolic Modeling

**"Model the metabolism of my microbial community"** -> Combine member genome-scale models into a community, then predict community growth, individual growth rates, and metabolite exchange (cross-feeding and competition) under a shared medium.
- Python: `micom.Community(taxonomy).cooperative_tradeoff()` (steady-state, abundance-weighted); SMETANA (cross-feeding scores); COMETS (dynamic)

## The governing principle: a community model inherits every member's errors, and the shared space is a modeling choice

Two things dominate whether a community prediction means anything:

- Reference-model quality propagates. A community model is only as good as its member reconstructions - a wrong biomass, a missing pathway, or an energy-generating cycle in one member distorts the whole community's exchange predictions. Curate members (systems-biology/model-curation) before combining, and confirm they share a namespace (BiGG vs ModelSEED); a namespace mismatch silently breaks metabolite sharing.
- How the shared space is modeled is the central design decision, and the classic trap is compartment pooling. Modeling a community as one giant "bag" with a single shared metabolite pool is fast but biologically wrong: it lets any member use any other member's INTERNAL metabolites directly, inventing cross-feeding that requires no secretion. The correct structure gives each member its own compartments and connects them only through a shared EXTRACELLULAR medium with explicit exchange. Pooling artifacts are a recurring reviewer catch. MICOM and SteadyCom implement the compartmentalized structure correctly; hand-merging models by prefix usually does not.

A further modeling fork: steady-state community FBA (SteadyCom, MICOM) assumes a stable coexistence with a common community growth rate, while dynamic simulation (COMETS, BacArena) resolves the time course and spatial structure but is expensive and parameter-hungry. Neither predicts the other's regime.

## Decision: which community method

| Goal | Tool | Approach / trade-off |
|------|------|----------------------|
| Metagenome-scale gut community, abundance-weighted, steady state | MICOM (Python) | community FBA with cooperative tradeoff (community vs individual growth); scales to many taxa from abundances |
| Cross-feeding / competition SCORES between members | SMETANA (CLI) | MRO (resource overlap = competition), MIP (interaction potential = cooperation), per-metabolite scores; pairs with CarveMe |
| Coexistence at a common community growth rate | SteadyCom | enforces one shared growth rate; elegant steady-state coexistence model |
| Time course / spatial dynamics, diffusion | COMETS / BacArena | dynamic (COMETS) or individual-based spatial (BacArena) FBA; realistic but expensive/parameter-hungry |

Do not model a community as one pooled "bag" model; use a tool that keeps members compartmentalized and connects them through a shared extracellular medium.

## Build and Simulate a Community with MICOM

**Goal:** Combine member models (weighted by their metagenomic abundance) and predict community and per-member growth under a medium.

**Approach:** Assemble a taxonomy table (one row per taxon with an `id`, a model `file`, and an `abundance`), build the `Community` (which compartmentalizes members correctly), and solve with cooperative tradeoff - which finds a community growth optimum while spreading growth across members rather than letting one taxon dominate. Reserve the `fraction` argument to trade community optimum against individual growth.

```python
from micom import Community
from micom.data import test_taxonomy

# taxonomy: columns id, file (per-taxon SBML), and abundance (from metagenomics). test_taxonomy()
# ships a ready E. coli example community.
taxonomy = test_taxonomy()

community = Community(taxonomy)            # builds the compartmentalized multi-species model
solution = community.cooperative_tradeoff(fraction=1.0)   # QP; needs HiGHS/CPLEX/Gurobi
print('community growth rate:', solution.growth_rate)
print(solution.members[['growth_rate']])  # per-taxon growth; NaN row is the shared medium
```

## Cross-Feeding and Competition (SMETANA)

```bash
# SMETANA (separate install) scores interactions between member models built by CarveMe:
#   pip install smetana   # then:
# smetana model1.xml model2.xml -o community --flavor bigg
# Outputs: MRO (metabolic resource overlap = competition for shared nutrients),
#          MIP (metabolic interaction potential = potential cooperation/cross-feeding),
#          and per-metabolite SMETANA scores (who feeds whom). A high MIP with low MRO
#          suggests cooperative cross-feeding; high MRO suggests competition.
```

## Dynamic and Spatial Simulation (COMETS)

```python
# For the time course rather than a steady state, COMETS (cometspy) runs dynamic FBA on a lattice
# with metabolite diffusion. Use when the QUESTION is temporal (succession, diauxie, spatial
# structure), not a coexistence steady state. It is far more expensive and needs kinetic parameters
# (uptake Vmax/Km, initial biomass, diffusion constants) that a steady-state model does not.
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Cross-feeding predicted that needs no secretion | compartment pooling (single shared internal pool) | use MICOM/SteadyCom (compartmentalized); connect members only via a shared extracellular medium |
| Members will not exchange metabolites | namespace mismatch (BiGG vs ModelSEED IDs) | reconcile member models via MetaNetX before combining |
| Community growth nonsensical | a member model is broken (bad biomass, energy cycle) | curate each member first; a bad member poisons the community |
| `cooperative_tradeoff` errors on solver | it is a QP and GLPK cannot solve it | use HiGHS (bundled), CPLEX, or Gurobi |
| One taxon takes all the growth | plain community-max FBA has alternate optima | use cooperative tradeoff (spreads growth) and set abundances from data |
| Dynamic run is impossibly slow | COMETS/BacArena are expensive and parameter-hungry | use a steady-state method unless the question is genuinely temporal/spatial |

## Related Skills

- systems-biology/metabolic-reconstruction - Build the member models (CarveMe pairs with SMETANA)
- systems-biology/model-curation - Curate members before combining; errors propagate to the community
- systems-biology/flux-balance-analysis - Single-organism FBA underlying each member
- metagenomics/abundance-estimation - Member abundances to weight the community
- metagenomics/functional-profiling - Community-level metabolic potential from metagenomes

## References

- Diener C, Gibbons SM, Resendis-Antonio O. 2020. MICOM: metagenome-scale modeling to infer metabolic interactions in the gut microbiota. *mSystems* 5(1):e00606-19.
- Zelezniak A, Andrejev S, Ponomarova O, et al. 2015. Metabolic dependencies drive species co-occurrence in diverse microbial communities. *PNAS* 112(20):6449-6454. (SMETANA)
- Chan SHJ, Simons MN, Maranas CD. 2017. SteadyCom: predicting microbial abundances while ensuring community stability. *PLoS Comput Biol* 13(5):e1005539.
- Zomorrodi AR, Maranas CD. 2012. OptCom: a multi-level optimization framework for the metabolic modeling and analysis of microbial communities. *PLoS Comput Biol* 8(2):e1002363.
- Harcombe WR, Riehl WJ, Dukovski I, et al. 2014. Metabolic resource allocation in individual microbes determines ecosystem interactions and spatial dynamics. *Cell Rep* 7(4):1104-1115. (COMETS)
- Dukovski I, Bajic D, Chacon JM, et al. 2021. A metabolic modeling platform for the computation of microbial ecosystems in time and space (COMETS). *Nat Protoc* 16(11):5030-5082.
- Bauer E, Zimmermann J, Baldini F, Thiele I, Kaleta C. 2017. BacArena: individual-based metabolic modeling of heterogeneous microbes in complex communities. *PLoS Comput Biol* 13(5):e1005544.
- Machado D, Andrejev S, Tramontano M, Patil KR. 2018. Fast automated reconstruction of genome-scale metabolic models for microbial species and communities. *Nucleic Acids Res* 46(15):7542-7553.

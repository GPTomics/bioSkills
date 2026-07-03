# Community Metabolic Modeling - Usage Guide

## Overview

Community metabolic modeling combines member genome-scale models into a multi-species model to predict community growth, individual growth rates, and metabolite exchange (cross-feeding and competition). Two facts govern whether the predictions mean anything. First, the community inherits every member's errors: a wrong biomass, missing pathway, or energy-generating cycle in one member distorts the whole community, and members must share a namespace or metabolite exchange silently fails. Second, how the shared space is modeled is the central design choice, and the classic mistake is compartment pooling - a single shared metabolite pool that lets any member use any other's internal metabolites, inventing cross-feeding that requires no secretion. Correct tools (MICOM, SteadyCom) keep members compartmentalized and connect them only through a shared extracellular medium. A separate fork is steady-state (SteadyCom/MICOM) versus dynamic/spatial (COMETS/BacArena) modeling.

## Prerequisites

```bash
pip install micom            # abundance-weighted steady-state community FBA; needs a QP solver (HiGHS bundled)
pip install smetana          # cross-feeding/competition scores (separate CLI; pairs with CarveMe)
# COMETS dynamic/spatial: cometspy + the COMETS engine (github.com/segrelab/comets)
```

Inputs: curated member models in SBML sharing one namespace, and member abundances (from metagenomics) for weighting. MICOM ships `test_taxonomy()` for a runnable example.

## Quick Start

Tell your AI agent:
- "Build a community model from these member SBMLs and predict community and per-member growth"
- "Score cross-feeding and competition between these two species with SMETANA"
- "Weight the community by my metagenomic abundances and run cooperative tradeoff"
- "Should I use MICOM steady-state or COMETS dynamic modeling for my question?"
- "Why is my hand-merged community inventing cross-feeding?"

## Example Prompts

### Community Growth
> "I have curated models for five gut species and their relative abundances. Build a MICOM community, run cooperative tradeoff, and report the community growth rate and each member's growth."

### Cross-Feeding
> "Use SMETANA on my two CarveMe models to compute MRO and MIP, and tell me whether these species mainly compete or cross-feed."

### Method Choice
> "I want to model succession over time in a co-culture with spatial structure - recommend steady-state MICOM versus dynamic COMETS and explain the parameter cost."

### Debugging Pooling
> "My community predicts a species growing on a metabolite no one secretes - check whether I built a pooled bag model and how to compartmentalize it properly."

## What the Agent Will Do

1. Confirm member models are curated and share a namespace; reconcile via MetaNetX if not.
2. Assemble a taxonomy table (id, model file, abundance) from metagenomics.
3. Build a compartmentalized community (MICOM) rather than a single-pool bag model.
4. Solve with cooperative tradeoff for community and per-member growth, or SMETANA for interaction scores.
5. For temporal/spatial questions, use COMETS/BacArena with the required kinetic parameters.
6. Flag pooling artifacts, namespace mismatches, and member-model defects that distort predictions.

## Tips

- Curate every member before combining; one broken member (bad biomass, energy cycle) poisons the whole community.
- Confirm all members share a namespace (BiGG vs ModelSEED); a mismatch means metabolites never connect and no cross-feeding appears.
- Never model a community as one shared metabolite pool ("bag" model) - it invents cross-feeding that needs no secretion. Use MICOM/SteadyCom, which compartmentalize members.
- Cooperative tradeoff (MICOM) spreads growth across members and needs a QP solver (HiGHS/CPLEX/Gurobi); plain community-max FBA has alternate optima and lets one taxon dominate.
- Weight members by measured abundances; an unweighted community answers a different question.
- SMETANA MRO measures competition (resource overlap), MIP measures cooperation potential; high MIP with low MRO suggests cross-feeding.
- Use steady-state methods unless the question is genuinely temporal or spatial; COMETS/BacArena need uptake kinetics, initial biomass, and diffusion constants a steady-state model does not.
- A cross-feeding prediction is a hypothesis; validate against measured metabolite exchange or co-culture data.

## Related Skills

- systems-biology/metabolic-reconstruction - Build member models (CarveMe pairs with SMETANA)
- systems-biology/model-curation - Curate members before combining; errors propagate
- systems-biology/flux-balance-analysis - Single-organism FBA underlying each member
- metagenomics/abundance-estimation - Member abundances to weight the community
- metagenomics/functional-profiling - Community metabolic potential from metagenomes

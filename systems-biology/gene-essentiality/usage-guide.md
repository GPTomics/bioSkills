# Gene Essentiality Analysis - Usage Guide

## Overview

In-silico gene essentiality predicts which genes a metabolic model cannot lose and still make biomass. The result is a conditional prediction about a specific model on a specific medium, never an intrinsic property of a gene. Three choices dominate the outcome: the growth cutoff that separates essential from viable (a policy you set, not a library default), the medium (rich vs minimal give different essential sets), and the deletion method (FBA re-optimizes the mutant, MOMA/ROOM model an immediate un-adapted mutant). Because genes reach reactions through gene-protein-reaction rules, a gene with an isozyme partner looks non-essential no matter how central its reaction - which is exactly why synthetic-lethal pairs exist. Predictions top out around 85-93% accuracy on E. coli metabolic genes; report MCC, not accuracy, and match the medium when validating against experimental screens.

## Prerequisites

```bash
pip install cobra pandas scikit-learn
# QP/MILP methods (MOMA, ROOM) benefit from an academic CPLEX or Gurobi; HiGHS also works.
```

Inputs: a genome-scale model (SBML/JSON or a built-in like `load_model('textbook')`) and, for validation, an experimental essential-gene list (Keio, Tn-seq, or CRISPR fitness) with the medium it was measured on.

## Quick Start

Tell your AI agent:
- "Find the essential genes in my model on glucose minimal medium"
- "Sweep the essentiality cutoff and tell me which calls are low-confidence"
- "Use MOMA instead of FBA because these are fresh transposon mutants"
- "Find synthetic-lethal gene pairs among the viable genes"
- "Compare my predicted essentials to the Keio set and report the MCC"

## Example Prompts

### Essential Genes
> "Run single-gene deletions on iML1515 in M9 glucose, call genes essential below 1% of wild-type growth, and also report how many calls change if I use 5% instead."

### Method Choice
> "These are freshly knocked-out CRISPR mutants scored in one growth cycle - use MOMA, not FBA, and explain how the essential set differs from the FBA call."

### Synthetic Lethality
> "Among genes that are individually non-essential, find synthetic-lethal pairs and flag any that could be combination drug targets."

### Condition-Specific
> "Compare essential genes aerobically vs anaerobically and separate the core-essential set from the condition-specific ones."

### Validation
> "Compare my model's predicted essential genes to the Keio collection on the same medium and report MCC, sensitivity, and the false-essential and false-non-essential genes."

## What the Agent Will Do

1. Load the model and confirm the medium matches the intended condition (or the experiment being compared to).
2. Take wild-type growth, then run `single_gene_deletion` (GPR-aware) and apply an explicit cutoff.
3. Sweep the cutoff to separate high-confidence essentials from sick-tail, cutoff-dependent calls.
4. For fresh mutants, use MOMA/ROOM against a wild-type reference instead of FBA.
5. For synthetic lethality, restrict to viable singles and run `double_gene_deletion`; use Fast-SL logic for higher-order sets.
6. Recompute essentiality per medium for condition contrasts.
7. Validate against experimental screens on the matched medium and report MCC.

## Tips

- Delete genes, not reactions - `single_gene_deletion` evaluates the GPR so complexes and isozymes are handled correctly.
- The `ids` column holds sets of gene-id strings; `list(s)[0]` gives the id (no `.id` attribute).
- State and sweep the cutoff. A result that flips between a 1% and 5% cutoff is a threshold artifact, not a finding.
- Match the medium to the experiment before comparing; a biosynthesis gene is essential on minimal medium and dispensable on rich.
- Use MOMA/ROOM for immediate/un-evolved mutants and FBA for evolved strains or when only hard lethality matters (the methods mostly agree on lethal-vs-viable).
- Cap the double-deletion gene list; the sweep is O(n^2). Use Fast-SL flux-support pruning for triples/quads.
- A central gene that comes back non-essential almost always has an isozyme - treat it as a synthetic-lethal candidate.
- False essentials usually mean the model is missing a bypass or isozyme (a curation problem); false non-essentials mean the real killer is non-metabolic (regulation, toxicity) and FBA cannot see it.
- Report MCC and sensitivity, not accuracy; essential genes are the minority class and accuracy is inflated by true negatives.
- Wrap deletion scripts in `if __name__ == '__main__':` (the functions spawn workers) or pass `processes=1`.

## Related Skills

- systems-biology/flux-balance-analysis - The FBA/medium/objective foundation for knockouts
- systems-biology/model-curation - Fix missing isozymes/bypasses that cause false essentials
- systems-biology/strain-design - Growth-coupling designs that build on knockout logic
- crispr-screens/hit-calling - Experimental essentiality screens for validation
- pathway-analysis/go-enrichment - Functional enrichment of essential-gene sets

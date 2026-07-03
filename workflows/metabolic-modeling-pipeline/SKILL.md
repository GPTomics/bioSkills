---
name: bio-workflows-metabolic-modeling-pipeline
description: End-to-end genome-scale metabolic modeling from genome sequence to flux predictions. Covers automated reconstruction with CarveMe, model validation with memote, FBA/FVA analysis, and gene essentiality prediction. Use when building metabolic models or predicting metabolic phenotypes from genomic data.
tool_type: mixed
primary_tool: cobrapy
goal_approach_exempt: true
workflow: true
depends_on:
  - systems-biology/metabolic-reconstruction
  - systems-biology/model-curation
  - systems-biology/flux-balance-analysis
  - systems-biology/gene-essentiality
  - systems-biology/context-specific-models
qc_checkpoints:
  - after_reconstruction: "Reactions 1000-2500, growth >0.01 on target media"
  - after_curation: "Memote score >50%, <5% orphan reactions"
  - after_fba: "Realistic growth rate, major pathways active"
  - after_essentiality: "Core essential genes match literature >70%"
---

## Version Compatibility

Reference examples tested with: COBRApy 0.29+, matplotlib 3.8+, numpy 1.26+, pandas 2.2+, seaborn 0.13+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Metabolic Modeling Pipeline

**"Build and analyze a metabolic model for my organism"** -> Orchestrate CarveMe reconstruction, memote quality scoring, gap-filling, FBA/FVA flux analysis, gene essentiality prediction, and context-specific model building from expression data.

Complete workflow for genome-scale metabolic modeling: from protein sequences to flux predictions and phenotype analysis.

## Workflow Overview

```
Protein FASTA (genome annotation)
        |
        v
[1. Reconstruction] --> CarveMe / gapseq / ModelSEED
        |
        v
[2. Model Curation] --> memote QC, gap-filling
        |
        | <---- Iterative refinement loop
        v
[3. FBA Analysis] --> Growth prediction, flux distribution
        |
        +-----------------------+
        |                       |
        v                       v
[4a. Gene Essentiality]   [4b. Context-Specific]
    Single/double KO       Tissue-specific models
        |                       |
        v                       v
Essential Gene List      Condition-Specific Fluxes
```

## Prerequisites

```bash
pip install cobra carveme memote escher pandas numpy matplotlib seaborn

conda install -c bioconda diamond
```

**Required data:**
- Protein FASTA file from genome annotation
- BiGG universal model (downloaded by CarveMe)

## Primary Path: Bacterial Model from Genome

### Step 1: Automated Reconstruction with CarveMe

```bash
# Basic reconstruction from a PROTEIN FASTA (CarveMe rejects raw/GenBank genomes)
carve genome.faa -o model_draft.xml

# Gram type / universe is a VALUE of -u/--universe, NOT a --gram-neg flag
carve genome.faa -o model_draft.xml -u gramneg

# Gap-fill for a specific medium (opt-in; the medium determines what gets added)
carve genome.faa -o model_draft.xml -u gramneg --gapfill M9
```

```python
import cobra

model = cobra.io.read_sbml_model('model_draft.xml')
print(f'Model: {model.id}')
print(f'Reactions: {len(model.reactions)}')
print(f'Metabolites: {len(model.metabolites)}')
print(f'Genes: {len(model.genes)}')

# Quick growth test
# Growth rate >0.01 h^-1 indicates viable model
solution = model.optimize()
print(f'Growth rate: {solution.objective_value:.4f} h^-1')
```

### Step 2: Model Validation with Memote

```bash
# Run the memote test suite (results stored as JSON when --filename is given)
memote run --filename model_result.json.gz model_draft.xml

# Generate the human-readable HTML snapshot report
memote report snapshot --filename model_report.html model_draft.xml
```

```python
# The memote SCORE measures consistency and annotation (well-formedness), NOT biological
# correctness -- a model can score high and mispredict every knockout. Read WHICH tests fail
# (stoichiometric consistency, mass/charge balance, energy-generating cycles) in the HTML report,
# and validate predictions separately (see systems-biology/model-curation). Programmatic access:
from memote.suite.api import test_model
code, result = test_model(model, results=True)   # result is a MemoteResult of the raw outcomes
```

### Step 3: Model Curation (Iterative)

```python
import cobra
from cobra.flux_analysis import gapfill

model = cobra.io.read_sbml_model('model_draft.xml')

# Check for common issues
def diagnose_model(model):
    issues = []

    # Dead-end metabolites (produced but not consumed, or vice versa)
    for met in model.metabolites:
        producing = [r for r in met.reactions if met in r.products]
        consuming = [r for r in met.reactions if met in r.reactants]
        if len(producing) > 0 and len(consuming) == 0:
            issues.append(f'Dead-end (not consumed): {met.id}')
        elif len(producing) == 0 and len(consuming) > 0:
            issues.append(f'Dead-end (not produced): {met.id}')

    # Blocked reactions
    fva = cobra.flux_analysis.flux_variability_analysis(model)
    blocked = fva[(fva['minimum'] == 0) & (fva['maximum'] == 0)]
    if len(blocked) > 0:
        issues.append(f'Blocked reactions: {len(blocked)}')

    return issues

issues = diagnose_model(model)
print(f'Found {len(issues)} issues')
for issue in issues[:10]:
    print(f'  {issue}')
```

```python
# Gap-filling for growth on specific media
from cobra.flux_analysis import gapfill

# Load universal reaction database for gap-filling
universal = cobra.io.read_sbml_model('universal_model.xml')

# Define target medium (e.g., glucose minimal)
target_medium = {
    'EX_glc__D_e': 10,  # Glucose uptake
    'EX_o2_e': 20,       # Oxygen
    'EX_nh4_e': 100,     # Ammonium
    'EX_pi_e': 100,      # Phosphate
    'EX_so4_e': 100,     # Sulfate
}

# Apply medium (model.exchanges yields Reaction objects, not id strings)
for rxn in model.exchanges:
    rxn.lower_bound = -target_medium[rxn.id] if rxn.id in target_medium else 0  # block other uptakes

# Gap-fill to enable growth
# Gap-filling adds minimal reactions from universal model to enable growth
gapfill_solution = gapfill(model, universal, demand_reactions=False)
print(f'Gap-fill added {len(gapfill_solution[0])} reactions')

# Gap-filled reactions are the LEAST-evidenced part of the model (added to force growth on this
# medium, not because homology supports them) -- flag them low-confidence, do not treat as validated.
for rxn in gapfill_solution[0]:
    model.add_reactions([rxn])
    print(f'  Added (low-confidence): {rxn.id} - {rxn.name}')

# Verify growth
solution = model.optimize()
print(f'Growth after gap-fill: {solution.objective_value:.4f} h^-1')
```

### Step 4: Flux Balance Analysis

```python
import cobra
import pandas as pd
import matplotlib.pyplot as plt

model = cobra.io.read_sbml_model('model_curated.xml')

# Basic FBA
solution = model.optimize()
print(f'Objective (growth): {solution.objective_value:.4f} h^-1')
print(f'Status: {solution.status}')

# Get active fluxes
fluxes = solution.fluxes
active_fluxes = fluxes[abs(fluxes) > 1e-6]
print(f'Active reactions: {len(active_fluxes)} / {len(model.reactions)}')

# Key exchange fluxes (uptake/secretion)
exchange_fluxes = fluxes[[r.id for r in model.exchanges]]
significant_exchanges = exchange_fluxes[abs(exchange_fluxes) > 0.1]
print('\nSignificant exchanges:')
print(significant_exchanges.sort_values())
```

```python
# Flux Variability Analysis (FVA)
from cobra.flux_analysis import flux_variability_analysis

# FVA identifies reaction flexibility
# Fraction 0.9 = allow 90% of optimal growth
fva = flux_variability_analysis(model, fraction_of_optimum=0.9)

# Identify rigid vs flexible reactions
fva['range'] = fva['maximum'] - fva['minimum']
rigid = fva[fva['range'] < 1e-6]
flexible = fva[fva['range'] > 1]

print(f'Rigid reactions (fixed flux): {len(rigid)}')
print(f'Flexible reactions: {len(flexible)}')

# Plot flux ranges for key pathways
glycolysis = ['PGI', 'PFK', 'FBA', 'TPI', 'GAPD', 'PGK', 'PGM', 'ENO', 'PYK']
glyc_fva = fva.loc[fva.index.isin(glycolysis)]

fig, ax = plt.subplots(figsize=(10, 6))
ax.barh(range(len(glyc_fva)), glyc_fva['maximum'] - glyc_fva['minimum'],
        left=glyc_fva['minimum'], alpha=0.7)
ax.set_yticks(range(len(glyc_fva)))
ax.set_yticklabels(glyc_fva.index)
ax.set_xlabel('Flux range (mmol/gDW/h)')
ax.set_title('Glycolysis Flux Variability')
plt.tight_layout()
plt.savefig('glycolysis_fva.pdf')
```

### Step 5a: Gene Essentiality Prediction

```python
from cobra.flux_analysis import single_gene_deletion, double_gene_deletion

# Single gene knockouts. Result columns: ids (a SET of gene-id strings), growth, status.
single_ko = single_gene_deletion(model)
single_ko['growth_ratio'] = single_ko['growth'] / solution.objective_value

# Essential genes: knockout drops growth below the cutoff (a policy, not a library default).
# Report and sweep the cutoff; match the medium to any experiment being compared. Essentiality
# is model- and medium-relative (see systems-biology/gene-essentiality).
essential = single_ko[single_ko['growth_ratio'] < 0.1]
print(f'Essential genes: {len(essential)} / {len(model.genes)} on this medium')

# ids elements are gene-id STRINGS (a set), so list(s)[0] gives the id -- there is no .id attribute.
essential_list = [list(s)[0] for s in essential['ids']]
with open('essential_genes.txt', 'w') as f:
    f.write('\n'.join(essential_list))
```

```python
# Double gene knockouts (synthetic lethality)
# WARNING: Computationally intensive for large models

# Focus on non-essential genes only (a synthetic lethal needs both singles viable)
non_essential = [g.id for g in model.genes if g.id not in essential_list]

# Run pairwise deletions (positional gene_list1/gene_list2; cap the O(n^2) sweep)
double_ko = double_gene_deletion(model, non_essential[:100], non_essential[:100])

# Synthetic lethality: neither single KO is lethal, but the double KO is
synthetic_lethal = double_ko[double_ko['growth'] < 0.01]
print(f'Synthetic lethal pairs: {len(synthetic_lethal)}')
```

### Step 5b: Context-Specific Models

Use a validated extraction method rather than ad-hoc pruning. COBRApy has no native GIMME/iMAT/INIT; the real Python options are troppo and corda, and the threshold/method choice dominates the result more than the data does. See systems-biology/context-specific-models for the method decision table and the threshold-sensitivity discipline.

```python
# corda is the most turnkey native-Python extraction method.
from corda import CORDA, reaction_confidence

# Translate expression into CORDA confidence classes (-1 absent, 0 unknown, 1 low, 2 med, 3 high)
# through the GPR, then build the context model.
gene_conf = {g.id: 2 for g in model.genes}                       # derive from expression quantiles
rxn_conf = {r.id: reaction_confidence(r, gene_conf) for r in model.reactions}   # pass the Reaction, not its GPR string
opt = CORDA(model, rxn_conf)
opt.build()
context_model = opt.cobra_model('tissue')
print(f'Context model: {len(context_model.reactions)} reactions')
# Rebuild at 2-3 thresholds and report which reactions are threshold-dependent (hypotheses).
```

## Visualization with Escher

```python
import escher

# Load model and solution
model = cobra.io.read_sbml_model('model_curated.xml')
solution = model.optimize()

# Create Escher map
builder = escher.Builder(
    map_name='e_coli_core.Core metabolism',
    model=model,
    reaction_data=solution.fluxes.to_dict()
)

builder.save_html('flux_map.html')
```

## Parameter Recommendations

| Step | Parameter | Value | Rationale |
|------|-----------|-------|-----------|
| CarveMe | --gapfill | M9 or LB | Match experimental media |
| Memote | score threshold | >50% | Minimum for usable model |
| FBA | solver | gurobi/cplex | Faster than glpk for large models |
| FVA | fraction_of_optimum | 0.9 | 90% allows realistic flexibility |
| Essentiality | growth threshold | 0.1 | Standard 10% of WT growth |
| Context | expression percentile | 25 | Balance specificity vs viability |

## Troubleshooting

| Issue | Likely Cause | Solution |
|-------|--------------|----------|
| No growth | Missing essential reactions | Gap-fill with universal model |
| Unrealistic growth rate | Unbounded uptake | Constrain medium properly |
| Many blocked reactions | Dead-end metabolites | Check metabolite connectivity |
| Low memote score | Missing GPR, mass balance | Run memote report for details |
| Essentiality mismatch | Missing isozymes | Add alternative pathways |

## Output Files

| File | Description |
|------|-------------|
| `model_draft.xml` | Initial reconstruction (SBML) |
| `model_curated.xml` | Gap-filled and validated model |
| `model_report.html` | Memote QC report |
| `essential_genes.txt` | Predicted essential genes |
| `fba_fluxes.tsv` | Optimal flux distribution |
| `fva_results.tsv` | Flux variability ranges |
| `flux_map.html` | Escher visualization |

## Extensions

Beyond the core genome-to-flux path, the model feeds two further analyses: build a multi-species community from several reconstructions (systems-biology/community-metabolic-modeling), or design growth-coupled knockouts to overproduce a target chemical (systems-biology/strain-design).

## Related Skills

- systems-biology/metabolic-reconstruction - CarveMe, gapseq details
- systems-biology/model-curation - Memote, gap-filling, energy-generating-cycle checks
- systems-biology/flux-balance-analysis - FBA, FVA, pFBA, sampling
- systems-biology/gene-essentiality - Single/double knockouts, MOMA/ROOM
- systems-biology/context-specific-models - Tissue-specific models (troppo/corda)
- systems-biology/community-metabolic-modeling - Multi-species community FBA (MICOM/SMETANA)
- systems-biology/strain-design - Growth-coupled knockout design (OptKnock/RobustKnock)

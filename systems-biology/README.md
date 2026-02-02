# systems-biology

## Overview

Constraint-based metabolic modeling including flux balance analysis, genome-scale model reconstruction, curation, and context-specific model building.

**Tool type:** mixed | **Primary tools:** cobrapy, CarveMe, memote, gapseq

## Skills

| Skill | Description |
|-------|-------------|
| flux-balance-analysis | FBA and FVA for metabolic flux prediction with COBRApy |
| metabolic-reconstruction | Build draft models from genomes with CarveMe, gapseq |
| model-curation | Validate and gap-fill models with memote |
| gene-essentiality | In silico gene knockouts and synthetic lethality |
| context-specific-models | Tissue-specific models with GIMME, iMAT algorithms |

## Example Prompts

- "Run FBA on the E. coli core model"
- "Predict growth rate on glucose minimal media"
- "Find essential genes in my metabolic model"
- "Build a metabolic model from this genome sequence"
- "Create a liver-specific model using GTEx expression data"
- "Identify synthetic lethal gene pairs"

## Requirements

```bash
pip install cobra escher
# CarveMe requires diamond and cplex/gurobi (or use glpk for open-source)
pip install carveme memote
```

## Related Skills

- **pathway-analysis** - Pathway enrichment context
- **metabolomics** - Integrate with metabolomics data
- **metagenomics** - Community metabolic models

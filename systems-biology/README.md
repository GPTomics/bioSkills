# systems-biology

## Overview

Constraint-based / genome-scale metabolic modeling: draft reconstruction from genomes, curation and gap-filling, FBA/FVA flux prediction, gene essentiality, context-specific model building, multi-species community modeling, and computational strain design.

**Tool type:** mixed | **Primary tools:** cobrapy, CarveMe, gapseq, memote, MICOM, StrainDesign

## Skills

| Skill | Description |
|-------|-------------|
| metabolic-reconstruction | Build draft genome-scale models from genomes with CarveMe / gapseq |
| model-curation | Validate, gap-fill, and standardize models with memote (consistency vs predictive validity) |
| flux-balance-analysis | Predict growth and flux distributions (FBA/FVA/pFBA/sampling) with COBRApy |
| gene-essentiality | In-silico single/double knockouts, synthetic lethality, FBA vs MOMA/ROOM |
| context-specific-models | Tissue/condition-specific models via GIMME/iMAT/CORDA (troppo, corda) |
| community-metabolic-modeling | Multi-species community FBA and cross-feeding with MICOM / SMETANA |
| strain-design | Growth-coupled knockout design with StrainDesign (OptKnock/RobustKnock/MCS) |

## Example Prompts

- "Build a metabolic model from this bacterial genome"
- "Curate my model and check it can't make ATP from nothing"
- "Run FBA and predict growth on glucose minimal media"
- "Find the essential genes in my model on M9"
- "Build a liver-specific model from my expression data"
- "Model my gut microbiome community and predict cross-feeding"
- "Design knockouts to overproduce succinate, coupled to growth"

## Requirements

```bash
pip install cobra memote carveme micom straindesign escher
# CarveMe needs diamond + an LP solver (CPLEX/Gurobi academic, or glpk/HiGHS open-source)
# gapseq is cloned from GitHub (github.com/jotech/gapseq); SMETANA/COMETS are separate installs
```

## Related Skills

- **metabolomics** - Integrate measured metabolite/flux data (isotope-tracing for 13C-MFA)
- **pathway-analysis** - Pathway/functional context for model genes
- **metagenomics** - Member abundances and functional potential for community models
- **genome-annotation** - Produce the annotated protein FASTA CarveMe/gapseq consume

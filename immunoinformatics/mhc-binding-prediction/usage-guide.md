# MHC Binding Prediction - Usage Guide

## Overview

Predict peptide-MHC class I binding affinity and natural-presentation likelihood to nominate candidate CD8 T-cell epitopes, with the methodological judgment to read the right score (BA vs EL), threshold correctly (%Rank not nM), and know where the predictions fail (rare alleles, low-expression neoantigens). For CD4/HLA class II, see the mhc-class-ii-prediction skill.

## Prerequisites

```bash
pip install mhcflurry
mhcflurry-downloads fetch
# NetMHCpan-4.1 and MixMHCpred are standalone academic downloads (not pip);
# the IEDB REST API wraps NetMHCpan if a local install is unavailable.
```

## Quick Start

Tell your AI agent what you want to do:
- "Predict class I presentation for these peptides with my patient's HLA genotype"
- "Scan this protein for 9-mer epitopes and rank by %Rank"
- "Which of these neoantigen peptides are presented, accounting for expression?"
- "Should I read the binding-affinity or the eluted-ligand score for this question?"

## Example Prompts

### Presentation Prediction

> "Score GILGFVFTL and NLVPMVATV against HLA-A*02:01, A*24:02, B*07:02 and report the best-presenting allele"

> "Classify these peptides as strong/weak/non-binders using %Rank, not raw nM"

### Protein Scanning

> "Tile this spike protein into 8-11mers and return windows under 2% Rank for common HLA-A alleles"

> "Find candidate CD8 epitopes in this antigen for a specific patient genotype"

### Score Choice and Caveats

> "Is this a binding-affinity question or a presentation question, and which tool/score fits?"

> "These neoantigens are lowly expressed; how do I avoid the EL abundance bias misranking them?"

> "This allele is rare and non-European; how much should I trust the prediction?"

## What the Agent Will Do

1. Confirm the patient HLA class I genotype and whether the question is "can it bind" (BA) or "is it presented" (EL/presentation)
2. Choose a predictor (MHCflurry for scripting; NetMHCpan-4.1 for broadest coverage; MixMHCpred for MS-grounded motifs)
3. Score peptides in a batched call and report affinity (nM), %Rank, and presentation score
4. Threshold on %Rank (strong <= 0.5%, weak <= 2.0%) for cross-allele comparability
5. Flag rare/extrapolated alleles and integrate expression for neoantigen ranking
6. Caveat that a predicted binder is a candidate, not a validated epitope

## Tips

- **BA vs EL** - read binding-affinity for "can it physically bind"; read eluted-ligand/presentation for "is it naturally presented" (the default question)
- **%Rank, not nM** - the 500 nM convention is allele-biased; always threshold on %Rank for multi-HLA work
- **Abundance bias** - EL/MS models over-rank peptides from abundant proteins; integrate expression when scoring low-expression neoantigens
- **Rare alleles** - pan-models extrapolate and never say "I don't know"; verify training support before trusting a number
- **Length** - 8-11mers, 9mers dominate; non-9mers rest on thinner training data
- **Binding != epitope** - presentation is one stage of a multi-stage funnel; do not report it as immunogenicity

## Related Skills

- immunoinformatics/mhc-class-ii-prediction - CD4/HLA class II binding (open groove, register, DQ pairing)
- immunoinformatics/neoantigen-prediction - applies class I binding to tumor mutations
- immunoinformatics/immunogenicity-scoring - the separate, weaker prediction of T-cell response
- clinical-databases/hla-typing - determine the patient genotype

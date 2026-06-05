# MHC Class II Prediction - Usage Guide

## Overview

Predict peptide-HLA class II (DR/DQ/DP) binding and presentation to nominate candidate CD4 T-cell epitopes, with the judgment to handle the open-groove register problem, the DQ/DP heterodimer pairing trap, and the much lower reliability of class II versus class I. For CD8/HLA class I, see the mhc-binding-prediction skill.

## Prerequisites

```bash
# NetMHCIIpan-4.3 and MixMHC2pred-2.0 are standalone academic downloads (DTU / Gfeller lab),
# not pip-installable. Request a license, install the binary, and add it to PATH.
# The IEDB MHC-II REST API wraps NetMHCIIpan if a local install is unavailable.
pip install pandas   # for parsing tool output
```

## Quick Start

Tell your AI agent what you want to do:
- "Predict CD4 epitopes in this antigen for DRB1*01:01 and DRB1*15:01"
- "Score these long peptides against the patient's DQ and DP heterodimers"
- "Find class II neoantigens for CD4 help in my tumor"
- "How should I format DQ alleles, and how much should I trust the DQ calls?"

## Example Prompts

### CD4 Epitope Prediction

> "Tile this antigen into 15-mers and return NetMHCIIpan strong binders (Rank <= 1%) for DRB1*04:01"

> "Predict class II presentation for these peptides with MixMHC2pred and NetMHCIIpan and compare"

### Heterodimer Handling

> "Build the correct NetMHCIIpan allele strings for this DQA1/DQB1 genotype using documented pairings"

> "Should I expand all DQA1 x DQB1 combinations, or only cis pairs?"

### Reliability and Caveats

> "These are DQ-restricted calls; how confident should I be relative to DR?"

> "The reported binding core shifts between runs - is this register ambiguity?"

## What the Agent Will Do

1. Confirm the patient class II genotype (both alpha and beta chains for DQ/DP)
2. Format alleles correctly per tool (NetMHCIIpan vs MixMHC2pred differ; DR vs DQ/DP differ)
3. Restrict DQ/DP to documented cis pairings rather than the full combinatorial cross product
4. Score long peptides (12-25mers), reporting the inferred 9-mer core
5. Apply class II %Rank cutoffs (strong <= 1%, weak <= 5%), not class I cutoffs
6. Flag DQ/DP calls as lower-confidence than DR and treat all class II calls as ranked hypotheses

## Tips

- **Class II is less reliable** - report calls as ranked hypotheses, not facts; DR > DP > DQ in trustworthiness
- **DQ/DP need both chains** - they are alpha/beta heterodimers; DR is single-chain
- **Pairing trap** - do not expand all DQA1 x DQB1 combinations; use documented cis pairings and state the assumption
- **Looser thresholds** - strong <= 1%, weak <= 5% %Rank; copying class I 0.5%/2.0% over-filters
- **Register ambiguity** - the open groove allows multiple binding frames; unstable cores are expected, corroborate across tools
- **Nomenclature differs** - NetMHCIIpan `DRB1_0101` / `HLA-DQA10501-DQB10201`; MixMHC2pred `DRB1_15_01` / `DPA1_02_01__DPB1_01_01`

## Related Skills

- immunoinformatics/mhc-binding-prediction - CD8/HLA class I binding (the solved regime)
- immunoinformatics/neoantigen-prediction - class II neoantigens for CD4 help
- immunoinformatics/immunogenicity-scoring - CD4 immunogenicity prediction
- clinical-databases/hla-typing - resolve DR/DQ/DP alleles for both chains

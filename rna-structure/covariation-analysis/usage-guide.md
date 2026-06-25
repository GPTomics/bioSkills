# Covariation Analysis - Usage Guide

## Overview
Test whether a proposed or predicted RNA secondary structure is supported by evolutionary covariation using R-scape. R-scape scores compensatory substitutions across an alignment against a phylogeny-aware null and estimates statistical power, so the result is a three-way verdict (supports / rejects / cannot infer), not a simple pass/fail. This is the gold-standard validation that found no covariation support for the proposed HOTAIR, Xist, and SRA lncRNA structures.

## Prerequisites
```bash
# R-scape (includes the CaCoFold structure predictor and R2R diagram output)
conda install -c bioconda rscape

# Input: a Stockholm alignment (deep and diverse), with a #=GC SS_cons line for -s.
# Build alignments with Infernal cmalign, RNAalifold, or from an Rfam SEED.
```

## Quick Start
Tell your AI agent what you want to do:
- "Is this conserved RNA structure actually supported by covariation?"
- "Run R-scape on my alignment and tell me if the structure is real"
- "Does my alignment even have the power to test structure?"
- "Predict a covariation-supported structure from my alignment"
- "Validate my SS_cons before I build a covariance model"

## Example Prompts

### Validate a Proposed Structure
> "I have a Stockholm alignment with a consensus structure for a putative ncRNA. Run R-scape and tell me whether the covariation is significant or just phylogenetic."

> "My lncRNA has a published secondary structure. Test whether it has evolutionary support."

### Power and Interpretation
> "R-scape found no significant pairs. Is that evidence against the structure, or does my alignment lack power?"

> "How many diverse homologs do I need before a covariation test is meaningful?"

### De Novo Structure
> "I have no trusted structure. Build a covariation-supported consensus with CaCoFold."

### Feed Downstream
> "Validate this SS_cons with R-scape before I build a covariance model from it."

## What the Agent Will Do
1. Confirm the alignment is deep and diverse enough for a covariation test to have power
2. Run R-scape with `-s` to test the given SS_cons, or `--cacofold` to build a structure de novo
3. Read the covarying-pair counts and per-pair power from the `.cov` and `.power` outputs
4. Return the three-way verdict (supports / rejects / cannot infer), not just a pair count
5. Pass a validated or CaCoFold structure to covariance-model building or thermodynamic folding

## Tips
- **Power makes the negative meaningful** - "0 significant pairs" is only evidence against a structure if the alignment had power to detect covariation. A low-power negative is "cannot infer", not "rejects"; report which one.
- **Depth and diversity** - Covariation needs sequences that vary at paired columns while preserving the pair. Few or near-identical sequences give low power; gather more diverse homologs.
- **Significance, not score** - A positive covariation score is not validation; require the E-value to clear the target (default 0.05) against the phylogenetic null.
- **Alignment quality** - Misaligned columns destroy real covariation and can manufacture spurious signal; validate the alignment first.
- **CaCoFold for de novo** - When there is no trusted structure, `--cacofold` builds one from significant covariation (it can include pseudoknots), giving a strong consensus to seed a covariance model or compare against a thermodynamic fold.
- **The lncRNA lesson** - A thermodynamically plausible structure is not established until covariation is statistically demonstrated; many proposed lncRNA structures fail this test.
- **Keep outputs tidy** - Pass `--outdir` so the `.cov`/`.power`/R2R files do not land in the working directory.

## Related Skills
- secondary-structure-prediction - Predict the structure whose conservation is then tested
- ncrna-search - Validate a custom CM's SS_cons here before building the covariance model
- structure-probing - Experimental evidence complementary to evolutionary covariation
- alignment/msa-statistics - Assess the alignment depth and diversity covariation needs
- phylogenetics/tree-io - The phylogeny underlying the covariation null

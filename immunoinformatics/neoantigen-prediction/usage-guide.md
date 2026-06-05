# Neoantigen Prediction - Usage Guide

## Overview

Identify tumor neoantigens from somatic variants with pVACtools, with the methodological discipline the field actually requires: binding prediction is the easy part, and the analysis must center the downstream attrition (clonality, HLA loss-of-heterozygosity, expression, proximal-variant phasing, agretopicity/foreignness quality, and validation tiers), because a binding-only pipeline has single-digit-percent positive predictive value.

## Prerequisites

```bash
# pVACtools (current 7.x; dependencies are heavy - use a dedicated env)
conda create -n pvactools python=3.11
conda activate pvactools
pip install pvactools vatools
pvacseq install_vep_plugin $VEP_PLUGINS   # installs Wildtype + Frameshift
# Ensembl VEP, a local IEDB install, and an HLA typer (OptiType/arcasHLA/HLA-HD) are also needed.
```

## Quick Start

Tell your AI agent what you want to do:
- "Find neoantigens from my VEP-annotated somatic VCF for this patient's HLA"
- "Rank these neoantigens by tumor-specific quality, not just IC50"
- "Drop candidates on HLA alleles the tumor lost (LOHHLA)"
- "Incorporate proximal germline variants by phasing before prediction"

## Example Prompts

### End-to-End Calling

> "Run VEP with the Wildtype and Frameshift plugins, annotate expression, then run pVACseq with my class I alleles"

> "Call fusion-derived neoantigens from my Arriba output with pVACfuse"

### Prioritization and Quality

> "Compute agretopicity (DAI) correctly as the WT/MT binding ratio and flag anchor-position mutations"

> "Rank candidates by cancer cell fraction and foreignness, and mark subclonal ones"

### Tumor-Specific Pitfalls

> "Which candidates are presented by an allele lost to HLA LOH?"

> "These are frameshift neoantigens - how should I weight them given the MS-validation gap?"

## What the Agent Will Do

1. Confirm the VCF is VEP-annotated with Wildtype + Frameshift plugins and a protein FASTA
2. Annotate expression and read counts (VAtools) so the filters are not silent pass-throughs
3. Phase proximal somatic/germline variants and supply the phased VCF
4. Run pVACseq with the patient HLA and sane VAF/coverage/expression filters
5. Run LOHHLA and drop candidates on lost alleles; estimate purity/CCF for clonality
6. Rank by quality (agretopicity, foreignness, expression, clonality) and re-tier in pVACview
7. Frame the output as a tier-1 hypothesis list for MS and functional validation

## Tips

- **VEP plugins** - Wildtype + Frameshift (NOT Downstream, which was dropped in pVACtools 2.0)
- **Binding is the easy part** - spend effort on clonality, LOH, expression, and quality features
- **HLA LOH** - run LOHHLA and remove candidates on deleted alleles; it fails silently if skipped
- **Clonality** - use cancer cell fraction (purity + copy number), not raw VAF
- **Phasing** - supply `--phased-proximal-variants-vcf` or the predicted peptides are ones the tumor never makes
- **Predicted != presented != immunogenic** - the in-silico list is the input to validation, not the answer
- **Frameshifts/fusions** - high-value, high-foreignness, but under-validated and rarely seen by MS

## Related Skills

- immunoinformatics/mhc-binding-prediction - the binding step (the solved, low-leverage part)
- immunoinformatics/mhc-class-ii-prediction - class II neoantigens for CD4 help
- immunoinformatics/immunogenicity-scoring - quality ranking of the candidate list
- clinical-databases/hla-typing - the genotype substrate
- variant-calling/variant-calling - upstream somatic calls
- workflows/neoantigen-pipeline - end-to-end orchestration

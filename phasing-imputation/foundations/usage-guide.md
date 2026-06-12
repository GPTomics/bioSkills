# Phasing and Imputation Foundations - Usage Guide

## Overview

This is the conceptual spine of the phasing-imputation category. It frames the decisions that come before any engine runs: that phasing and imputation are one Li-Stephens copying HMM, that imputation's honest output is a dosage with a self-estimated quality rather than a hard genotype, that the pipeline stages are ordered and each fails silently when misaligned, and that the upstream choice of genotyping strategy (dense array vs low-coverage WGS) is an ascertainment decision, not just an accuracy one. It routes the per-stage mechanics to the four sibling skills (reference-panels, haplotype-phasing, genotype-imputation, imputation-qc) and carves the boundary to read-backed phasing, which is a physically different signal owned by long-read-sequencing.

Use this skill to plan a pipeline, choose a genotyping strategy, or diagnose imputation that runs clean but produces wrong or non-replicating results.

## Prerequisites

- bcftools (`conda install -c bioconda bcftools`) for the pre-flight readiness checks.
- A QC'd, biallelic input VCF on a known genome build (GRCh37 or GRCh38). The build and strand must be reconciled with the reference panel before phasing.
- Conceptual prerequisites worth holding in mind:
  - Imputation outputs dosages (expected alt-allele count in [0,2]) with a per-variant quality (INFO/R2/DR2); these must be filtered, and the dosage (not a hard call) is what downstream regression should use.
  - The pipeline order is fixed: QC -> align build and strand to the panel -> phase -> impute per chromosome -> filter by INFO/R2 plus a MAF floor -> dosage GWAS. Misordering fails silently.
  - The reference panel is the prior: imputation can only impute variants the panel contains, and accuracy degrades when panel ancestry does not match the target.
  - The downstream engines (Beagle, SHAPEIT5, Minimac4, IMPUTE5, GLIMPSE2) are Li-Stephens variants; the differences are engineering, not different theories.

## Quick Start

Tell your AI agent what you want to do:
- "Plan a phasing and imputation pipeline for my array cohort and tell me the stage order"
- "Should I use a SNP array plus imputation or low-coverage WGS for my admixed cohort?"
- "Check whether my VCF is ready to phase and impute"
- "My imputation ran fine but the GWAS hits do not replicate - what stage failed?"
- "Explain why I should carry dosages instead of hard-called genotypes into association testing"

## Example Prompts

### Strategy selection
> "I have a 20,000-sample admixed Latino cohort and a fixed budget. Walk me through whether a SNP array plus imputation against TOPMed or low-coverage (1x) WGS plus GLIMPSE2 is the better genotyping strategy, and explain the ascertainment trade-off."

### Pipeline sequencing
> "I have an Illumina array VCF on GRCh37. Lay out the full ordered pipeline from QC to a dosage-based GWAS, naming which sibling skill owns each stage and where a silent failure could occur."

### Diagnosing a silent failure
> "My Beagle imputation completed without errors, but accuracy is near zero across whole chromosomes and the GWAS lambda is off. Diagnose the likely build/strand/ordering problem."

### The dosage concept
> "Explain why an imputed genotype is a dosage and a self-estimated R2, not a measurement, and what that means for how I should run the association test and filter variants."

## What the Agent Will Do

1. Frame the request against the ordered pipeline and identify which stage(s) it touches.
2. For a strategy question, weigh array-plus-imputation vs low-coverage-WGS-plus-imputation by population ancestry, target variant frequency (common vs rare), and budget, then recommend a default.
3. Run a pre-flight readiness check (genome build, chr naming, biallelic fraction, strand-ambiguous SNP burden, MAF spectrum) before any phasing.
4. Route each stage to the owning sibling skill: panel choice to reference-panels, phasing to haplotype-phasing, imputation to genotype-imputation, quality filtering to imputation-qc.
5. Enforce the discipline that dosages (not hard calls) move downstream and that phase-dependent claims need a switch-error rate against an independent truth set.
6. Hand the filtered dosages to population-genetics/association-testing for the GWAS test.

## Tips

- Treat a confident-looking INFO/R2 in an under-represented ancestry with suspicion: the metric grades the model's own posterior and cannot see that the panel lacked the target's haplotypes.
- Record the panel name, genome build, and release date with every result; the same code against a different panel release imputes a different set of sites.
- Reach for read-backed (long-read/Hi-C) phasing when the goal is to phase a private variant or resolve a clinical compound heterozygote; statistical phasing cannot place a variant no one else carries.
- Pre-phase the cohort once and re-impute against newer panels without re-phasing, but invest in the phasing step, because its switch errors are baked into every later imputation.
- A genetic (recombination) map is in centimorgans on a specific build; a wrong-build or flat map degrades phasing silently. Match the map to the data build.

## Related Skills

- reference-panels - Select and prepare the panel; ancestry-match, build, strand alignment, server access
- haplotype-phasing - Statistical phasing engines and the switch-error metric
- genotype-imputation - Imputation engines, dosage fields, servers, low-coverage WGS
- imputation-qc - INFO/R2/DR2 metrics, MAF-stratified accuracy, dosage filtering
- long-read-sequencing/haplotype-phasing - Read-backed / molecular single-sample phasing
- variant-calling/vcf-basics - Source and basic handling of the input VCF
- population-genetics/association-testing - GWAS test on the imputed dosages
- population-genetics/population-structure - PCA to establish target ancestry for panel selection
- workflows/gwas-pipeline - End-to-end QC -> phase -> impute -> associate -> visualize

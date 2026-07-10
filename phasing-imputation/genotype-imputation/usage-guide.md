# Genotype Imputation - Usage Guide

## Overview

Impute untyped genotypes against a reference panel - for array data (Beagle, Minimac4, IMPUTE5) or low-coverage WGS from genotype likelihoods (GLIMPSE2, QUILT2, STITCH). The idea that shapes downstream use is that an imputed genotype is a posterior, not a measurement: the deliverable is a dosage (expected alt-allele count in [0,2]) plus a self-estimated quality (Beagle DR2, Minimac R2, IMPUTE INFO) that estimates r2 from the posterior spread without ever seeing the truth. Downstream GWAS regresses on dosages, not hard calls. This skill covers engine choice, the dosage/quality output fields, the phasing prerequisite, per-chromosome chunking, the Michigan/TOPMed servers, and low-coverage WGS as the modern array replacement. It routes phasing, panel preparation, and quality filtering to sibling skills.

## Prerequisites

- An imputation engine: Beagle (a Java jar), Minimac4, IMPUTE5, or GLIMPSE2 for low-coverage WGS (`conda install -c bioconda` for several). bcftools for FORMAT-field extraction and genotype-likelihood generation.
- A reference panel in the engine's format (bref3 for Beagle, msav for Minimac4, imp5 for IMPUTE5), prepared and build/strand-aligned (route to reference-panels), and a build-matched genetic map.
- For Minimac4/IMPUTE5, a PHASED target (route to haplotype-phasing); Beagle and GLIMPSE2 phase internally.
- Conceptual prerequisites and big notes:
  - Output dosages (DS), not hard calls, go to association testing; request DS/HDS/GP explicitly.
  - The quality field (DR2/R2/INFO) is an estimate of imputation quality, not a measured accuracy, and cannot detect panel-ancestry mismatch.
  - Impute all samples together; separate case/control imputation manufactures false associations.
  - Run per chromosome; HRC and TOPMed are server-only.
  - Low-coverage WGS imputation consumes genotype likelihoods (PL/GL), not hard calls.

## Quick Start

Tell your AI agent what you want to do:
- "Impute my phased array VCF against TOPMed and give me dosages with the imputation quality"
- "Impute my low-coverage (1x) WGS BAMs with GLIMPSE2 against a reference panel"
- "Which imputation engine should I use for a very large reference panel run locally?"
- "Extract dosages and DR2 from my imputed VCF for filtering before GWAS"
- "Prepare my data for the Michigan Imputation Server"

## Example Prompts

### Array imputation
> "I have a phased Illumina array VCF for 5,000 European samples on GRCh38. Impute it against an appropriate panel per chromosome with Beagle, output dosages and the per-variant quality, and tell me which field carries the dosage."

### Low-coverage WGS
> "I have 0.5x whole-genome sequencing on 2,000 samples and a TOPMed-like panel. Run the GLIMPSE2 chunk/split/phase/ligate workflow from genotype likelihoods and explain why low-coverage WGS beats an array for my admixed cohort."

### Engine and server choice
> "My cohort needs HRC. Should I run Minimac4 locally or use the Michigan Imputation Server, and what are the access and reproducibility implications?"

### Dosage handling
> "Explain why I should carry dosages (DS) rather than hard-called genotypes into my GWAS, and how DR2/R2/INFO differ across Beagle, Minimac, and IMPUTE5."

## What the Agent Will Do

1. Confirm the input is phased (for Minimac4/IMPUTE5) or note that Beagle/GLIMPSE2 phase internally.
2. Choose the engine by data type (array vs low-coverage WGS), panel access, and scale.
3. Impute per chromosome against the prepared, build/strand-aligned panel, requesting DS/HDS/GP.
4. For low-coverage WGS, generate or read genotype likelihoods and run the GLIMPSE2 chunk/split-reference/phase/ligate pipeline.
5. For controlled-access panels (HRC, TOPMed), route to the imputation server and prepare per-chromosome uploads.
6. Hand dosages and the quality field to imputation-qc for filtering, then to population-genetics/association-testing.

## Tips

- Carry dosages (DS) into the regression; hard-calling discards exactly the uncertainty imputation exists to quantify.
- Treat DR2/R2/INFO as an estimate of imputation quality, not a validated accuracy; the only true accuracy comes from masking known genotypes (route to imputation-qc).
- Impute cases and controls together; batch-differential imputation quality is a classic source of non-replicating GWAS hits.
- Use the correct Minimac4 syntax (positional args, `minimac4 panel.msav target.vcf.gz`); the old `--refHaps`/`--haps` style is Minimac3.
- For under-represented ancestries or rare-variant work, consider low-coverage WGS plus GLIMPSE2 instead of an array; it removes the array's ascertainment bias.
- Request the FORMAT fields up front (`-f GT,DS,HDS,GP` for Minimac4; `gp=true ap=true` for Beagle), or the dosage will be missing downstream.

## Related Skills

- haplotype-phasing - Pre-phasing the target (required by Minimac4 and IMPUTE5)
- reference-panels - Select and prepare the panel and align build/strand
- imputation-qc - Filter by DR2/R2/INFO and MAF; the metric is an estimate, not truth
- variant-calling/vcf-basics - Genotype likelihoods (PL/GL) for low-coverage imputation
- variant-calling/variant-normalization - Split multiallelics before imputation
- population-genetics/association-testing - GWAS test on the imputed dosages
- clinical-databases/polygenic-risk - Polygenic scores from imputed dosages
- workflows/gwas-pipeline - End-to-end QC -> phase -> impute -> associate

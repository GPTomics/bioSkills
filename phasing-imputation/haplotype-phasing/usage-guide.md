# Haplotype Phasing - Usage Guide

## Overview

Estimate haplotype phase from population linkage disequilibrium - turning unphased genotypes (0/1) into phased haplotypes (0|1) using SHAPEIT5, SHAPEIT4, Eagle2, or Beagle. The idea that governs use is that statistical phase is an inference, not a measurement: it works well for common variants in LD with their neighbors and degrades steeply for rare variants, so the deliverable is a switch-error rate stratified by minor allele count, not a single genome-wide number, and rare-variant cis/trans calls (compound heterozygotes) need biobank-scale phasing or orthogonal trio/read-backed evidence. This skill owns the statistical paradigm and carves the boundary to read-backed phasing, which is a physically different signal owned by long-read-sequencing.

## Prerequisites

- A phasing engine: SHAPEIT5 (`conda install -c bioconda shapeit5`), Eagle2, or Beagle (a Java jar). bcftools for normalization and inspection.
- A QC'd, biallelic VCF/BCF and a genetic (recombination) map matching the data's genome build.
- A reference panel if reference-based phasing a small cohort (route to reference-panels for selection).
- Conceptual prerequisites and big notes:
  - Only heterozygous sites carry phase ambiguity; all error is measured there.
  - SHAPEIT5 is a suite (phase_common, phase_rare, ligate, switch), not one command; this changed from SHAPEIT4.
  - The genetic map must match the genome build; a mismatched map degrades phasing silently.
  - chrX male non-PAR is haploid and must be coded so; split multiallelics before phasing.
  - Reference-based phasing wins for small cohorts; within-cohort phasing wins past tens of thousands of samples.

## Quick Start

Tell your AI agent what you want to do:
- "Phase my array VCF with Eagle2 against an ancestry-matched reference panel"
- "Phase my biobank WGS data including rare variants with SHAPEIT5"
- "Phase my genotypes with Beagle as input to imputation"
- "Benchmark my phasing switch-error rate against a trio, stratified by allele frequency"
- "Should I phase against a reference panel or within my cohort?"

## Example Prompts

### Pre-phasing for imputation
> "I have a 5,000-sample European array VCF on GRCh38. Phase it as input to imputation against a matched panel, per chromosome, and tell me whether reference-based or within-cohort phasing is appropriate at this sample size."

### Rare-variant phasing
> "I have 50,000 WGS samples and need to call compound heterozygotes. Run the SHAPEIT5 common-scaffold-then-rare pipeline and explain why rare-variant phasing needs this design and how to report accuracy by minor allele count."

### Benchmarking
> "I have parent-parent-child trios in my cohort. Use them to compute my phaser's switch-error rate stratified by minor allele count and explain what switch versus Hamming error tells me."

### Boundary
> "I have long-read sequencing on a single sample and want to phase a private variant. Should I run SHAPEIT, and how does read-backed phasing fit in?"

## What the Agent Will Do

1. Confirm the genome build and align the genetic map to it; normalize the VCF to biallelic.
2. Choose reference-based vs within-cohort phasing by cohort size and panel availability, and the engine (SHAPEIT5 for rare-variant biobank work, Eagle2/Beagle for common-variant pre-phasing).
3. For rare-variant phasing, run phase_common -> ligate -> phase_rare with overlapping chunks.
4. Handle chrX by coding male non-PAR haploid and splitting the PARs.
5. Benchmark with the `switch` tool against trios where available, reporting SER stratified by minor allele count.
6. Hand the phased haplotypes to genotype-imputation, or flag when a phase-dependent claim needs trio/read-backed confirmation.

## Tips

- Do not trust a single genome-wide switch-error rate for a rare-variant call; phasing quality is a steep function of minor allele count.
- A switch error is invisible to genotype QC because it changes which haplotype an allele sits on, not the genotype; report a rate against an independent truth set.
- Make chunked phasing regions overlap so the ligate step can resolve phase across the seam; abutting chunks guarantee a switch.
- Combine, do not rank, the phase sources: seed SHAPEIT with read-backed phase from long reads when phasing a private variant, since statistical phasing cannot place a variant no one else carries.
- Match the genetic map and reference panel to the data build; a GRCh37 map on GRCh38 data corrupts phase with no error.

## Related Skills

- reference-panels - Select the ancestry-matched panel that reference-based phasing copies from
- genotype-imputation - Imputation consumes the phased haplotypes (pre-phasing)
- imputation-qc - Switch-error benchmarking alongside imputation quality QC
- long-read-sequencing/haplotype-phasing - Read-backed / molecular single-sample phasing
- variant-calling/variant-normalization - Split multiallelics and left-align before phasing
- causal-genomics/fine-mapping - Phased haplotypes feed haplotype-level fine-mapping
- clinical-databases/hla-typing - HLA typing is a high-stakes consumer of long-range phase
- workflows/gwas-pipeline - End-to-end QC -> phase -> impute -> associate

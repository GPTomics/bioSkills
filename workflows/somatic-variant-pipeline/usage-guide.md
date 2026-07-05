# Somatic Variant Pipeline Usage Guide

## Overview

An end-to-end tumor-normal somatic pipeline chains a somatic caller (Mutect2 or Strelka2) with the four somatic-specific removers - panel of normals, gnomAD germline-resource prior, cross-sample contamination estimate, and the FFPE/oxoG orientation-bias model - then somatic SV/CNV, then interpretation. The governing reality: there is no universal somatic truth set (the DREAM/Alioto benchmark showed low cross-pipeline concordance), somatic variants sit at continuous sub-1.0 VAF set by purity/ploidy/clonality, and interpretation uses the AMP/ASCO/CAP tier system plus oncogenicity - never germline ACMG.

## Prerequisites

```bash
conda install -c bioconda gatk4 strelka manta bcftools ensembl-vep
```

Resources: a genome reference (`.fa` + `.dict` + `.fai`), af-only-gnomAD VCF (germline prior), a small common biallelic-SNP VCF such as `small_exac_common_3.vcf.gz` (contamination), and a 40+ normal panel of normals matched to the assay.

## Quick Start

Tell your AI agent what you want to do:
- "Call somatic mutations from my tumor-normal BAM pair"
- "Run Mutect2 with contamination and orientation-bias filtering"
- "Build a panel of normals from my normal samples for the same assay"
- "Call somatic variants tumor-only and warn me about germline leakage"
- "Run Mutect2 and Strelka2 and take the consensus"
- "How do purity and copy number change how I read the VAF?"
- "Which tier and oncogenicity framework applies to my somatic variants?"

## Example Prompts

### Calling
> "Call somatic SNVs and indels from tumor.bam vs normal.bam with Mutect2, including the full contamination and orientation-bias filtering chain."

> "Set up a panel of normals from 40 normal BAMs sequenced on the same platform, then use it in Mutect2."

### Tumor-only and reproducibility
> "I only have the tumor (archival FFPE, no matched normal) - call somatic variants tumor-only and tell me what false positives to expect."

> "Run Mutect2, Strelka2, and MuSE and give me the 2-of-3 consensus somatic callset."

### Purity and interpretation
> "My tumor is low purity - why are expected driver mutations missing, and how does purity/ploidy change VAF?"

> "Annotate my somatic VCF and route it to AMP/ASCO/CAP tiers and oncogenicity, not germline ACMG."

## What the Agent Will Do

1. Confirm inputs (tumor + matched normal BAMs, reference, af-only-gnomAD, common-SNP resource, PoN) or set up tumor-only with caveats.
2. Run Mutect2 (tumor+normal in one command) with the PoN and germline-resource priors, emitting `--f1r2-tar-gz`.
3. Learn the orientation-bias model (FFPE/oxoG) and estimate contamination on common biallelic SNPs.
4. Apply FilterMutectCalls with the contamination table, tumor segmentation, and orientation priors; extract PASS.
5. Optionally add a Strelka2 arm (with Manta candidateSmallIndels) and build a multi-caller consensus.
6. Reason about VAF against purity/ploidy and copy number; add somatic SV/CNV and TMB/MSI/signatures.
7. Normalize, annotate (VEP/Funcotator), and route to tier/oncogenicity interpretation.

## Tips

- Always use a matched normal when available - tumor-only has a materially higher false-positive rate and cannot cleanly separate somatic from germline.
- Build the PoN from 40+ normals on the SAME platform/chemistry; a mismatched PoN imports the wrong artifact profile.
- GetPileupSummaries needs a COMMON biallelic-SNP resource, not the af-only-gnomAD used for the germline prior.
- Always run the orientation-bias model for FFPE/archival input - it removes the C>T/G>A (FFPE) and C>A/G>T (oxoG) low-VAF artifacts.
- Low-purity tumors need more depth to hit the same sensitivity; get purity/ploidy from copy-number analysis and correct VAF to cancer-cell fraction before calling a variant subclonal.
- Normalize every caller's VCF (`bcftools norm -f ref.fa -m-`) before consensus intersection or annotation, or shared calls silently drop.
- Never apply germline ACMG to a tumor variant - use AMP/ASCO/CAP tiers (tumor-type-specific) and oncogenicity.

## Related Skills

- variant-calling/gatk-variant-calling - Mutect2 mechanism and germline context
- variant-calling/filtering-best-practices - FilterMutectCalls internals and normalization
- variant-calling/variant-annotation - VEP/Funcotator/SnpEff, transcript choice, COSMIC
- variant-calling/clinical-interpretation - AMP/ASCO/CAP tiers and oncogenicity
- variant-calling/structural-variant-calling - Somatic SV detection with Manta/GRIDSS
- copy-number/cnvkit-analysis - Somatic CNV and purity/ploidy

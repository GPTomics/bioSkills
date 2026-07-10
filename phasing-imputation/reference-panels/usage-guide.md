# Reference Panels - Usage Guide

## Overview

Select and prepare the reference panel that statistical phasing and imputation copy haplotypes from. The decision that dominates results is ancestry match: imputation can only impute variants the panel contains and copies haplotypes from samples that resemble the target, so a panel whose ancestry does not match the study population imputes poorly no matter how large it is. This skill covers panel choice (1000 Genomes, HRC, TOPMed, HGDP+1kGP, CAAPA), genome-build and chromosome-naming reconciliation, the strand/allele harmonization gate that silently corrupts results when skipped, and conversion to the engine-specific panel format. It does not run the phasing or imputation engines (those are sibling skills) and routes classical HLA-allele panels out to clinical-databases.

## Prerequisites

- bcftools (`conda install -c bioconda bcftools`) with the `fixref` and `fill-tags` plugins; PLINK for the allele-frequency file the harmonization check consumes.
- Will Rayner's `HRC-1000G-check-bim.pl` (or the bgen/VCF variant) for the field-standard strand/allele check, plus the panel's site/legend file.
- A reference FASTA matching the PANEL's genome build, and a QC'd, normalized study VCF.
- Conceptual prerequisites and big notes:
  - A panel is data with a dated release and a fixed build (GRCh37 or GRCh38). Record the exact version and build; "1000 Genomes" alone is unreproducible.
  - HRC is SNP-only (no indels) with a MAF floor near 5e-4; TOPMed is the most diverse but is server-only and never downloadable.
  - Match the panel build to the study build and avoid liftover; if forced to lift, do it once with a strand-aware method and re-run the harmonization check.
  - Normalize (split multiallelics, left-align) before aligning alleles to the panel.

## Quick Start

Tell your AI agent what you want to do:
- "Which reference panel should I use for my admixed Latino cohort on GRCh38?"
- "Align my GRCh37 array data to the HRC panel strand and alleles before imputation"
- "Convert my reference VCF to Beagle bref3 and Minimac msav formats"
- "My study data is GRCh37 but I want a GRCh38 panel - how should I handle the build?"
- "Check whether my data can use TOPMed given our data-residency rules"

## Example Prompts

### Panel selection
> "My cohort is an admixed US population sequenced on GRCh38. Which reference panel should I impute against and why, and note any access or governance requirements."

### Strand and build harmonization
> "I have an Illumina array VCF on GRCh37. Walk me through aligning strand and alleles to the panel, handling palindromic SNPs, and reading the allele-frequency concordance plot before I upload to an imputation server."

### Format conversion
> "I have a custom reference panel as a phased VCF. Convert it to bref3 for Beagle and msav for Minimac4, per chromosome, and tell me what each engine requires."

### Governance constraint
> "Our IRB forbids participant genotypes leaving the country. Which diverse panels can I use locally, and what accuracy trade-off does that imply versus TOPMed?"

## What the Agent Will Do

1. Establish the study genome build and the target ancestry (routing to population-genetics/population-structure for the PCA if needed).
2. Recommend a panel by ancestry match, build, variant class (SNP vs indel), and data-governance constraints, naming the exact version.
3. Normalize the study VCF and reconcile chromosome naming to the panel convention.
4. Run the strand/allele harmonization check, execute the generated fix script, drop unresolvable palindromic SNPs, and inspect the allele-frequency concordance plot against the ancestry-matched sub-panel.
5. Convert the panel to the engine-specific format (msav, bref3, or imp5), pairing it with a build-matched genetic map.
6. Hand the harmonized study data and prepared panel to haplotype-phasing and genotype-imputation.

## Tips

- When someone says "I used the biggest panel," the real question is how many panel samples share the target ancestry, because that subset did the work.
- Treat "all palindromic SNPs kept" as evidence the strand was never checked; the check drops A/T and C/G SNPs above MAF 0.4 deliberately.
- Run the harmonization check and then actually execute its generated fix script; the check diagnoses, the script fixes.
- A high INFO/R2 from an ancestry-mismatched panel can be confidently wrong; the metric is panel-relative and cannot detect that the panel lacked the target's haplotypes.
- If the data cannot be uploaded to a US server, the practical choice narrows to downloadable panels (1000G, HGDP+1kGP); HGDP+1kGP is the strong diverse-and-local default.
- Obtain both the site-only legend (for the check) and the full haplotype panel (for imputation); they are different files.

## Related Skills

- haplotype-phasing - The phasing engine that consumes the panel
- genotype-imputation - Impute untyped variants once the panel is prepared
- imputation-qc - INFO/R2 quality metrics, which cannot detect ancestry mismatch
- variant-calling/variant-normalization - Split multiallelics and left-align before harmonization
- population-genetics/population-structure - PCA to establish target ancestry
- clinical-databases/hla-typing - Classical HLA-allele imputation with a dedicated panel
- workflows/gwas-pipeline - End-to-end QC -> phase -> impute -> associate

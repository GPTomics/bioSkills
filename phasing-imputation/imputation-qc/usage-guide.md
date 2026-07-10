# Imputation QC - Usage Guide

## Overview

Assess and filter phasing and imputation output. The core idea is that the routine quality score (Beagle DR2, Minimac R2, IMPUTE/GLIMPSE INFO) is the model grading its own posterior, not a measured accuracy: it estimates r2 from the dosage variance without ever seeing the truth, it is confounded with allele frequency, and a flat INFO >= 0.3 cutoff is therefore a hidden rare-variant filter. True accuracy comes only from masking known genotypes and computing dosage-r2 binned by MAF; concordance lies for rare variants. The second trap is the differential-imputation confound: imputing cases and controls separately manufactures false GWAS hits that every per-group QC metric passes. This skill covers the metrics, MAF-stratified filtering, masked-truth accuracy, phasing switch-error QC, and dosage-based downstream usage.

## Prerequisites

- bcftools (`conda install -c bioconda bcftools`) for filtering and field extraction; cyvcf2, pandas, and numpy for the Python QC summaries; PLINK2 for downstream dosage handling.
- An imputed VCF with the engine's quality field (DR2, R2, or INFO) and AF; a truth/array genotype set for masked-accuracy validation; trios for phasing switch-error QC.
- Conceptual prerequisites and big notes:
  - The quality score is an estimate of r2 from the posterior, not a validation; it cannot detect panel-ancestry mismatch.
  - The field name is engine-specific (DR2 Beagle, R2 Minimac, INFO GLIMPSE/IMPUTE), and the numbers are not comparable across engines.
  - Pair any quality cutoff with a MAF floor, and prefer MAF-stratified filtering.
  - Impute cases and controls together; the differential-imputation confound is structural, not fixable by a filter.
  - Carry dosages, not hard calls, into association.

## Quick Start

Tell your AI agent what you want to do:
- "Filter my imputed VCF by DR2 and MAF before GWAS"
- "Report imputation quality stratified by minor allele frequency"
- "Validate my imputation accuracy by masking typed genotypes"
- "Check my phasing switch-error rate against trios"
- "Why do my imputed case-control hits not replicate?"

## Example Prompts

### Quality filtering
> "I have a Beagle-imputed VCF. Filter it for GWAS using the correct quality field plus a MAF floor, and explain why a single flat INFO cutoff acts as a hidden rare-variant filter."

### Accuracy validation
> "I have array genotypes and an imputed VCF. Mask the typed sites, compute the true dosage-r2 binned by minor allele frequency, and tell me why I should not report concordance for rare variants."

### Diagnosing inflation
> "My imputed case-control GWAS has genome-wide-significant hits that do not replicate and the QC looks clean. Could separate imputation of cases and controls be the cause, and how do I confirm it?"

### Phasing QC
> "I have parent-parent-child trios. Compute my phaser's switch-error rate stratified by minor allele count and explain switch versus Hamming error."

## What the Agent Will Do

1. Detect the engine's quality field (DR2, R2, or INFO) and extract it with AF.
2. Filter on the correct field with a MAF floor, MAF-stratified where rare variants matter, and state the thresholds.
3. Where a truth/array set exists, mask typed genotypes, re-impute, and compute dosage-r2 binned by MAF as the gold-standard accuracy curve.
4. Check the Rsq-versus-EmpRsq gap and the allele-frequency concordance plot for strand/ancestry problems (routing strand fixes to reference-panels).
5. Benchmark phasing switch-error against trios where available, stratified by minor allele count.
6. Confirm imputation was done jointly across batches, and hand filtered dosages to population-genetics/association-testing.

## Tips

- Use the engine's actual field: DR2 for Beagle, R2 for Minimac, INFO for GLIMPSE/IMPUTE; a filter on the wrong name silently passes everything.
- `bcftools +fill-tags` computes AF/MAF/HWE but cannot produce an imputation quality score; that number comes only from the imputer.
- A negative Minimac EmpR is a strand/allele flip; a large Rsq-versus-EmpRsq gap flags panel, strand, or ancestry mismatch.
- Concordance is disqualified for rare variants because a do-nothing imputer scores ~98% at MAF 1%; use masked dosage-r2 by MAF bin.
- The absence of a QC red flag never clears the differential-imputation confound; impute together and verify quality does not differ by batch.

## Related Skills

- genotype-imputation - Produces the dosages and the quality field this skill filters
- haplotype-phasing - Switch-error benchmarking of the phasing that precedes imputation
- reference-panels - Panel ancestry mismatch, which the quality metric cannot detect
- variant-calling/vcf-statistics - Generic VCF INFO/FORMAT field parsing
- population-genetics/association-testing - Consumes the filtered dosages
- long-read-sequencing/haplotype-phasing - Read-backed phasing switch QC
- workflows/gwas-pipeline - End-to-end QC -> phase -> impute -> associate

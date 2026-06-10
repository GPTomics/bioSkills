# Array QC and Filtering

## Overview

This skill decides which probes and samples on an Illumina Infinium methylation array (450K, EPIC v1, or EPIC v2) can be trusted, and how to merge data safely across array versions. It takes the corrected beta/M matrix produced upstream and removes probes whose values are noise, genotype, or cross-hybridization artifacts, collapses EPICv2 replicate probes to one value per locus, and runs the sample-level identity checks (sex prediction, rs-SNP fingerprinting) that catch the single most common array failure: a sample swap. It also diagnoses whether Sentrix chip and array position are confounded with the biological group.

The guiding idea is that a raw array beta is uninterpretable until it is detection-masked and probe-filtered. A failed, cross-reactive, or SNP-overlapping probe returns a confident-looking beta that is meaningless, and no downstream model can recover from a mislabeled sample or a chip-confounded batch. This is the measurement's integrity layer.

## Prerequisites

Install the R/Bioconductor stack:

```r
BiocManager::install(c('minfi', 'sesame', 'ChAMP',
                       'IlluminaHumanMethylation450kanno.ilmn12.hg19',
                       'IlluminaHumanMethylationEPICanno.ilm10b4.hg19'))
remotes::install_github('markgene/maxprobes')   # cross-reactive probe lists (not on Bioconductor)
```

Conceptual prerequisites:
- A corrected `GenomicRatioSet`/`RGChannelSet` (minfi) or beta matrix from array-preprocessing. This skill QC's a matrix it does not produce; detection-p and bead-count masking happen upstream.
- The array version must be known. Cross-reactive lists, SNP annotation, and the genome build are all array-version-specific.
- EPICv2 is hg38-native and requires sesame (or third-party manifest/annotation packages); 450K and EPIC v1 are hg19.
- Annotation packages are large and array-specific; install the one matching the data.

## Quick Start

Tell your AI agent what you want to do:
- "Filter SNP-overlapping and cross-reactive probes from my 450K data"
- "Collapse the EPICv2 replicate probes before I test per-CpG"
- "Check every sample's predicted sex against my sample sheet"
- "Cluster samples on the rs-SNP probes to find duplicates and swaps"
- "Merge my 450K and EPICv2 cohorts safely"

## Example Prompts

### Probe filtering
> "I have a corrected GenomicRatioSet from a 450K EWAS. Drop probes that overlap a SNP at the CpG or SBE site, remove the Chen 2013 cross-reactive list, drop chrX/chrY probes, and report how many probes survive each step."

### EPICv2 replicate collapse
> "My EPICv2 beta matrix has duplicated cg IDs with design suffixes. Collapse the replicate probes to one value per locus, keeping the lowest-detection-p replicate, before I run differential testing."

### Sample identity QC
> "Before I analyze this EPIC cohort, predict each sample's sex and flag any that disagree with the sample sheet, then cluster the samples on the rs-SNP genotyping probes to find unintended duplicates or swaps."

### Cross-version harmonization
> "I want to combine a 450K cohort with an EPICv2 cohort. Collapse the EPICv2 replicates, intersect the shared probes, and liftover so the coordinates are on the same build before merging."

### Batch diagnosis
> "Run an SVD diagnostic to tell me whether Sentrix chip or array position loads on the top variance components of my beta matrix, so I know if the batch is confounded with my groups."

## What the Agent Will Do

1. Confirms the array version (450K / EPIC v1 / EPIC v2) and that the input matrix has already been detection-masked and corrected upstream; routes back to array-preprocessing if not.
2. For EPICv2, collapses replicate probes to one value per cg core ID with `betasCollapseToPfx` before any per-locus step.
3. Drops SNP-overlapping probes with `dropLociWithSnps(snps=c('CpG','SBE'), maf=)`, choosing the maf cutoff for the question.
4. Removes the array-version-matched cross-reactive list via maxprobes.
5. Decides sex-chromosome handling: drop chrX/chrY for autosomal EWAS, or keep and sex-stratify when sex is the phenotype.
6. Runs sample-level identity QC: `getSex()` prediction versus the sample sheet, and rs-SNP fingerprint clustering with `getSnpBeta()`.
7. Diagnoses chip/position batch association (ChAMP `champ.SVD()` or equivalent) and hands explicit correction to ewas-design.
8. For cross-version merges, intersects probe IDs after collapse and lifts coordinates with `mLiftOver`.
9. Reports probe attrition at each step for the methods section.

## Tips

- Always start from raw IDATs or a properly masked matrix from array-preprocessing, never from an unmasked beta matrix where the failed-probe information has been discarded.
- The cross-reactive probe list is array-version-specific: Chen 2013 for 450K, Pidsley 2016 for EPIC. Using the wrong list drops the wrong probes.
- Over-filtering is a real cost. Dropping every flagged probe can remove clock or signature CpGs; treat the maf cutoff and cross-reactive list as calibrated, not automatic, and confirm target CpGs survive.
- A `getSex()` mismatch or two samples clustering together on rs-SNP betas is a swap until proven otherwise. Resolve it before any analysis.
- If chip or array position is confounded with the biological group, no filtering rescues it. That is a design problem for ewas-design; the fix is randomization at sample layout time.
- EPICv2 is hg38 and 450K/EPIC are hg19. Track the build per array version and liftover before merging coordinates or applying an older-array clock.

## Related Skills

- array-preprocessing - Produces the corrected beta/M matrix being QC'd and filtered
- ewas-design - Chip/position batch correction and study design
- cell-type-deconvolution - Cohort composition QC
- differential-cpg-testing - Downstream per-CpG testing on the filtered matrix
- workflows/methylation-pipeline - End-to-end pipeline

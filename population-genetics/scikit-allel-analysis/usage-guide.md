# scikit-allel Analysis - Usage Guide

## Overview

scikit-allel is the in-memory Python toolkit for population genetics, built on three array types (GenotypeArray, HaplotypeArray, AlleleCountsArray), but its statistics are mostly ratios and densities with one denominator bug in two faces: per-base diversity (pi, theta, Dxy) silently divides by total span instead of accessible bases unless an `is_accessible=` callability mask is passed, and a genome-wide FST must be a ratio-of-sums (`sum(a)/(sum(a)+sum(b)+sum(c))`), not a mean of per-SNP FST - scikit-allel deliberately returns the variance components so the aggregation is done correctly. Both bugs run silently and return plausible numbers, so correctness is a choice the analyst makes.

## Prerequisites

- scikit-allel installed: `pip install scikit-allel` (plus `zarr` and `dask[array]` for out-of-core work)
- Input genotypes as VCF/BCF, with phased haplotypes if selection scans (iHS, XP-EHH, nSL) are planned
- A callable-loci / accessibility mask (a per-base boolean from coverage and mappability) for any per-base diversity statistic
- A confident ancestral allele (outgroup alignment or an AA field) only if an unfolded SFS is needed; otherwise the folded SFS is the safe default
- Conceptual prerequisites and big notes:
  - Per-base pi/theta/Dxy need `is_accessible=`; without it the denominator is total span, not callable bp, deflating values and distorting the landscape.
  - A genome-wide FST is `sum(a)/(sum(a)+sum(b)+sum(c))` (or `sum(num)/sum(den)`), never `mean(per_snp_fst)`; `average_hudson_fst` does the aggregation plus a block-jackknife SE.
  - `to_n_alt()` defaults to `fill=0`, which imputes every missing call as homozygous reference and biases PCA/LD toward the reference allele.
  - `sfs()` is unfolded and expects DERIVED counts, not ALT counts; use `sfs_folded(ac)` when polarization is uncertain.
  - iHS/XP-EHH/nSL require phased haplotypes and return UNSTANDARDIZED scores; iHS/nSL standardize binned by derived-allele frequency (standardize_by_allele_count), XP-EHH genome-wide (standardize).
  - `read_vcf` is eager and OOMs on whole genomes; convert once with `vcf_to_zarr` and use `GenotypeDaskArray` for scale.
  - scikit-allel is in maintenance mode; sgkit is the named successor but is not yet at feature parity, so scikit-allel remains the pragmatic choice for established stat workflows.

## Quick Start

Tell your AI agent what you want to do:
- "Calculate nucleotide diversity from my VCF using an accessibility mask"
- "Compute genome-wide Hudson FST between two populations with a jackknife standard error"
- "Run a Patterson f3 admixture test to check whether my target population is admixed"
- "LD-prune my genotypes and run a Patterson-scaled PCA"
- "Scan phased haplotypes for selection with iHS and standardize the scores"
- "Convert my large VCF to Zarr and compute statistics out-of-core"

## Example Prompts

### Diversity and the SFS
> "Compute pi and Watterson's theta in 100kb windows using a callable-loci mask, and report the accessible base count per window."

> "Build a folded site frequency spectrum because I don't trust my ancestral-allele calls."

### Population differentiation
> "Estimate genome-wide FST between population A and population B with a block-jackknife standard error, using the Hudson estimator."

> "Give me a per-window FST landscape across the chromosome and flag the top outlier windows."

### Admixture
> "Run a Patterson f3 test with my target population as C to see if it is admixed between A and B."

> "Compute an ABBA-BABA D-statistic with a jackknife z-score for my four populations."

### Structure and selection
> "LD-prune my genotypes over three rounds and run a randomized Patterson PCA, then plot PC1 vs PC2."

> "Compute iHS across phased chromosome 2 and standardize the scores binned by derived-allele frequency."

### Large data
> "Convert my whole-genome VCF to Zarr and count alleles out-of-core with dask."

## What the Agent Will Do

1. Load the VCF region into a GenotypeArray with only the needed fields, and derive allele counts (missing-aware) and per-population counts.
2. For per-base diversity, require and apply an accessibility mask, and report the accessible base count per window.
3. For differentiation, aggregate FST as a ratio-of-sums (or use `average_hudson_fst`) with a block-jackknife SE, never a mean of per-SNP FST.
4. For admixture, run f3 with the test population first or a D-statistic, and interpret the jackknife z-score against the conventional thresholds.
5. For structure, build a missing-free 012 matrix, LD-prune iteratively, mask known inversions, and run Patterson-scaled PCA.
6. For selection, reshape phased haplotypes, run the scan, and standardize binned by derived-allele frequency.
7. At biobank scale, convert to Zarr once and use dask-backed arrays for out-of-core counting and filtering.

## Tips

- The single most common silent bug is omitting `is_accessible=`; the mask must come from coverage/mappability (a callable-loci BED), never from the VCF's variant positions.
- Prefer the `average_*` FST and f-statistic functions: they do the correct ratio-of-sums and a block-jackknife SE in one call.
- Size the jackknife block (`blen`) to exceed the LD decay length; too-small blocks correlate and produce anticonservative standard errors and false significance.
- Put the test population FIRST in f3 (`average_patterson_f3(acc, aca, acb, blen)`); a significantly negative f3 (z < ~-3) is the admixture signal.
- LD-prune over ~3 rounds and mask large polymorphic inversions by position; a single inversion survives pruning and will dominate a PC.
- Use `sfs_folded` whenever the ancestral allele is uncertain; the plain `sfs()` is unfolded and silently treats ALT as derived.
- For greenfield biobank-scale infrastructure, evaluate sgkit, but expect scikit-allel to remain necessary where sgkit lacks feature parity.

## Related Skills

- selection-statistics - selection-scan design, standardization, and demography-aware outlier interpretation (this skill owns the array mechanics)
- population-structure - PCA and ADMIXTURE via PLINK2/FlashPCA2
- linkage-disequilibrium - LD pruning and clumping
- plink-basics - PLINK-format QC before array-based analysis
- variant-calling/vcf-basics - VCF generation and manipulation before loading
- phasing-imputation/haplotype-phasing - phased haplotypes required for iHS/XP-EHH/nSL

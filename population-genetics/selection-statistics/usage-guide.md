# Selection Statistics - Usage Guide

## Overview

Selection statistics scan genomes for natural selection, but no single statistic separates selection from demography at one locus: a bottleneck mimics a sweep in the site-frequency spectrum, recent expansion gives genome-wide negative Tajima's D, and population structure gives positive D that mimics balancing selection. The honest deliverable is therefore empirical genome-wide outliers, multiple orthogonal signals (haplotype plus SFS plus diversity plus differentiation), and explicit demographic-model calibration, never an absolute cutoff. iHS and nSL detect ongoing/incomplete sweeps and collapse to zero at fixation, while XP-EHH and Rsb catch fixed sweeps by borrowing a second population, so they are complementary rather than redundant.

## Prerequisites

- scikit-allel and numpy: `pip install scikit-allel numpy`
- selscan 2.0 and the companion `norm` binary for EHH scans: `conda install -c bioconda selscan`
- SweepFinder2 (and optionally SweeD, OmegaPlus) for SFS-CLR scans
- Phased haplotypes plus a genetic map for iHS/XP-EHH; an ancestral-allele assignment for any derived-allele test
- Conceptual prerequisites and big notes:
  - No single statistic is a selection test; rank windows and flag the genome-wide empirical tail (top 1%), then intersect orthogonal signals.
  - iHS/nSL are standardized within derived-allele-frequency bins; XP-EHH/iHH12 get a plain genome-wide z-score - conflating the two is a real bug.
  - Derived-allele tests (Fay & Wu H, Zeng E, unfolded SFS) need substitution-model polarization with two or more outgroups, not a single chimp allele.
  - Background selection mimics FST, PBS, and CLR peaks in low-recombination regions; deconfound CLR with a SweepFinder2 B-value map.
  - Haplotype statistics require phased data and degrade with switch error; use selscan 2.0 `--unphased` if phasing is unreliable.

## Quick Start

Tell your AI agent what you want to do:
- "Scan for selection using genome-wide outliers, not absolute cutoffs"
- "Compute iHS for ongoing sweeps and standardize within frequency bins"
- "Run XP-EHH to find completed sweeps in one of two populations"
- "Compute Hudson FST and PBS for population-differentiation outliers"
- "Find soft sweeps with Garud's H12"
- "Run a SweepFinder2 CLR scan with a background-selection map"
- "Intersect outliers across iHS, FST, and reduced diversity"

## Example Prompts

### SFS Neutrality Tests
> "Compute windowed Tajima's D and nucleotide diversity, but interpret them relative to the genome-wide distribution rather than an absolute cutoff."

> "Run a Fay & Wu H scan and tell me how to polarize ancestral alleles correctly so CpG sites don't fake a sweep."

### Haplotype Tests
> "Compute iHS scores for ongoing selection and standardize them within derived-allele-frequency bins."

> "Run XP-EHH between my two populations to find completed sweeps, and apply the correct genome-wide standardization."

> "My data is poorly phased - can I still run a haplotype sweep scan?"

### Differentiation
> "Calculate Hudson FST between two populations with unequal sample sizes and report the ratio-of-averages."

> "Compute PBS across three populations to find which branch the differentiation is on."

> "An FST peak sits in a low-recombination region - is it selection or background selection?"

### Multi-Statistic and Polygenic
> "Intersect outliers from iHS, FST, and depressed diversity to get credible sweep candidates."

> "I want to test for polygenic adaptation on a trait - what effect sizes should I use?"

## What the Agent Will Do

1. Clarify the question (localized recent sweep vs ancient/recurrent coding selection vs polygenic trait shift) and route coding selection to comparative-genomics and polygenic shifts to Qx.
2. Assess data: phase status, genetic map availability, and whether ancestral alleles are polarized with an adequate outgroup model.
3. Pick complementary statistics by sweep completeness (iHS/nSL for incomplete, XP-EHH/Rsb for fixed, H12 for soft, CLR to localize, FST/PBS for differentiation).
4. Standardize correctly: derived-allele-frequency bins for iHS/nSL, genome-wide z-score for XP-EHH/iHH12.
5. Flag empirical genome-wide outliers (top 1%) rather than absolute cutoffs, and intersect orthogonal signals.
6. Deconfound background selection (B-value map for CLR, recombination control for FST) and calibrate against a demographic-model simulation where possible.
7. Report candidate regions with coordinates and the evidence each statistic contributes.

## Tips

- A null iHS is not evidence of no selection; a completed sweep erases the within-population contrast iHS needs - switch to XP-EHH/Rsb.
- Report the windowed proportion of |iHS|>2 SNPs, not single extreme scores, which are noisy.
- Filter rare variants before any FST/PBS outlier scan and use the Hudson estimator for unequal sample sizes (Bhatia 2013).
- nSL needs no genetic map and is more robust to recombination-rate variation; it is a good default when no reliable map exists.
- Conservation scores (phyloP, GERP, phastCons) measure purifying constraint, the opposite sign of a recent positive sweep - do not read them as "under selection."
- Any polygenic Qx result is only as trustworthy as the stratification control in the GWAS supplying its effect sizes; prefer within-family/sibling estimates (the height-selection retraction).
- One significant statistic is a hypothesis; the deliverable is the intersection of orthogonal signals plus demographic-model calibration.

## Related Skills

- scikit-allel-analysis - genotype/haplotype array loading and allele-count basics
- population-structure - PCA and ADMIXTURE for population assignment before FST/PBS
- linkage-disequilibrium - EHH and LD mechanics underlying haplotype statistics
- phasing-imputation/haplotype-phasing - phased haplotypes that all EHH statistics require
- comparative-genomics/positive-selection - dN/dS and McDonald-Kreitman for recurrent coding selection
- comparative-genomics/introgression-detection - archaic/admixture signals that confound differentiation outliers

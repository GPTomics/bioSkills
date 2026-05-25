# Conservation Genetics Usage Guide

## Overview

Assesses genetic health of populations for conservation management with effective-population-size estimation across multiple time horizons (LD-based contemporary Ne via NeEstimator V2 option-file API with physical-linkage correction via SNeP for genomic data, recent Ne trajectory over ~200 generations via GONE/GONE2 with mandatory populated cM map column, deep demographic history via PSMC for single high-coverage genomes or Stairway Plot 2 / dadi / fastsimcoal2 from SFS), F-statistics (hierfstat Weir-Cockerham with bootstrap CIs), runs of homozygosity binned by length class to date inbreeding events (bcftools roh HMM and detectRUNS), genetic-load decomposition into realized (homozygous deleterious) and masked (heterozygous deleterious) components per Bertorelle 2022, the modern 100/1000 Ne rule (Frankham 2014; revised from 50/500), taxon-specific Ne/Nc ratio with explicit caveats for high-fecundity species (Hauser & Carvalho 2008), purging dynamics per Hedrick & Garcia-Dorado 2016, tree-sequence-recorded forward simulations via SLiM with pyslim+tskit, and the Sukumaran-Knowles caveat that multispecies-coalescent methods are inappropriate for management-unit definition.

## Prerequisites

- R with hierfstat (`install.packages('hierfstat')`)
- adegenet for VCF/genepop/PLINK ingestion (`install.packages('adegenet')`)
- detectRUNS for ROH from PLINK (`install.packages('detectRUNS')`)
- poppr for private alleles (`install.packages('poppr')`)
- NeEstimator V2 (Java; https://github.com/bunop/NeEstimator2.X)
- GONE2 (`remotes::install_github('esrud/GONE2')` or standalone)
- SNeP (https://sourceforge.net/projects/snepnetrends/) for physical-linkage-corrected LDNe
- Stairway Plot 2 (Java; https://github.com/xiaoming-liu/stairway-plot-v2)
- PSMC (https://github.com/lh3/psmc)
- bcftools >= 1.19 for HMM-based ROH detection
- SLiM 4 (https://messerlab.org/slim/) with pyslim/tskit for forward simulations

## Quick Start

Tell your AI agent what you want to do:

- "Estimate contemporary Ne using NeEstimator V2 (option-file API, not CLI flags)"
- "Estimate recent Ne trajectory with GONE2; verify the cM map column is populated first"
- "Run F_ROH analysis and bin ROH by length class to date inbreeding events"
- "Apply SNeP physical-linkage correction to LDNe on my RAD-seq data"
- "Decompose genetic load into realized vs masked components (Bertorelle 2022)"
- "Reconstruct demographic history with PSMC OR Stairway Plot 2"
- "Simulate forward with SLiM 4 using tree-sequence recording for speed"

## Example Prompts

### Contemporary Ne with Method Choice

> "Estimate contemporary Ne from my single-sample SNP data via NeEstimator V2. Configure the option file properly (line order matters; defaults can silently fail). For RAD-seq, apply SNeP physical-linkage correction or thin SNPs to >= 1 cM apart."

### Recent Ne Trajectory

> "Run GONE2 on my PLINK BED/BIM/FAM data to detect recent bottlenecks over the last 200 generations. CRITICAL: verify the BIM cM column is populated (not zero); without this GONE2 will silently produce nonsense."

### Inbreeding via ROH

> "I have PLINK ped/map files from a captive breeding program. Detect ROH using `bcftools roh -G30 --AF-tag AF` (HMM-based, density-independent), bin by length class: > 16 Mb = parents/grandparents inbreeding, 1-4 Mb = ~10-20 generations back, < 1 Mb = deep background."

### Genetic Load Decomposition

> "Annotate my VCF with SnpEff. Decompose genetic load into realized (homozygous high-impact) and masked (heterozygous high-impact) components per Bertorelle 2022. Reference Hedrick & Garcia-Dorado 2016 for purging vs drift interpretation."

### Demographic History

> "Run PSMC on my whole-genome BAM (>= 20x coverage). Use mutation rate 1.4e-8 and generation time 5 years. For populations without WGS, switch to Stairway Plot 2 from the SFS."

> "Run fastsimcoal2 with >= 50 independent optimization replicates on my multi-population SFS. Report best-likelihood replicate AND the spread; single-replicate inference is unreliable due to local optima."

### Forward Simulation

> "Simulate a non-Wright-Fisher demographic scenario with selection in SLiM 4 using tree-sequence recording. Recapitate with msprime for the neutral burn-in and add neutral mutations with msprime.sim_mutations a posteriori."

### Ne/Nc Caveat

> "Convert my estimated Ne to census size N. Do NOT default to Ne/Nc = 0.1; cite Hauser & Carvalho 2008 for taxon-specific ratios (Ne/Nc 2-6 orders of magnitude smaller than census, i.e., 10^-2 to 10^-6, documented for marine fish with sweepstakes recruitment; 0.1-0.3 for long-lived mammals)."

## What the Agent Will Do

1. Identify which time horizon the Ne question requires and select the matching estimator (contemporary LD, recent trajectory GONE2, deep PSMC/Stairway Plot 2, multi-population dadi/fastsimcoal2)
2. For NeEstimator V2: construct the `.ne2` option file with strict line ordering; confirm methods actually ran via the INFO output line
3. For GONE2: verify the PLINK BIM cM column is populated before running; if not, populate from species genetic map or Hi-C-derived recombination map
4. For LDNe on genomic data: apply SNeP physical-linkage correction OR thin SNPs to >= 1 cM apart; do not naively run LDNe on linked SNPs
5. For ROH: use `bcftools roh -G30 --AF-tag AF` (HMM-based, density-independent) for VCF input or `detectRUNS` for PLINK; bin by length class to date inbreeding events
6. For genetic load: decompose realized vs masked per Bertorelle 2022; cite Hedrick & Garcia-Dorado 2016 for purging vs drift dynamics
7. For deep history: PSMC for single high-coverage diploid; Stairway Plot 2 from SFS without parametric model; dadi/fastsimcoal2 with >= 50 optimization replicates and report convergence
8. For forward simulation: SLiM 4 with `treeSeqOutput()`; recapitate with msprime for the neutral burn-in; add neutral mutations a posteriori via msprime.sim_mutations
9. Cite the modern 100/1000 Ne rule (Frankham 2014), NOT the obsolete 50/500
10. Cite Hauser & Carvalho 2008 when discussing Ne/Nc; do NOT default to 0.1
11. Cite the Sukumaran-Knowles 2017 oversplitting caveat when defining management units; do not use BPP/BFD* species-delimitation methods for ESU/MU work

## Tips

- The 50/500 rule has been revised to 100/1000 by Frankham et al. 2014 (Ne >= 100 for short-term, Ne >= 1000 for long-term adaptive potential)
- Ne/Nc varies dramatically by life history; the 0.1 default is for terrestrial vertebrates and is wrong for marine fish (Hauser & Carvalho 2008 documented Ne/Nc spanning 2-6 orders of magnitude smaller than census, i.e., 10^-2 to 10^-6, in sweepstakes-recruitment species) and other high-fecundity taxa
- NeEstimator V2 is option-file driven; CLI flags do not work; the `.ne2` line order matters and silent failures are common
- GONE2 silently produces nonsense if the BIM/MAP cM column is zero (PLINK default uses physical position); populate first
- For RAD-seq or WGS, naive LDNe overestimates drift due to physical linkage; use SNeP (Barbato 2015) or thin to >= 1 cM apart
- F_ROH is more informative than F_IS for individual inbreeding; bin ROH by length to date the inbreeding (> 16 Mb = recent, < 1 Mb = ancient)
- bcftools roh uses an HMM and works at lower SNP density than detectRUNS' window approach; prefer bcftools roh for RAD-seq with sparse coverage
- PSMC requires >= 20x WGS coverage; for lower coverage or RAD-seq, switch to Stairway Plot 2 from the SFS
- dadi and fastsimcoal2 likelihood surfaces have local optima; >= 50 independent replicates are required to report reliable demographic inference; single-replicate inference is methodologically incomplete
- SLiM 4 with `treeSeqOutput()` is 5-100x faster than mutation tracking; pyslim and tskit add neutral mutations after-the-fact
- Genetic load is NOT a single number; decompose into realized (homozygous deleterious) and masked (heterozygous deleterious) per Bertorelle 2022; cite Hedrick & Garcia-Dorado 2016 for purging-vs-drift dynamics
- Robinson 2018 (Channel Island foxes) is a key empirical example of purging; Kyriazis 2021 provides simulation-based load assessment
- Frankham 2015 meta-analysis shows genetic rescue benefits outweigh outbreeding-depression risk in most documented cases; the asymmetry favors action
- Do NOT use BPP / BFD* / multispecies-coalescent methods to define management units; Sukumaran-Knowles 2017 showed these oversplit, mistaking population structure for species
- Always pair Ne estimates with census-size context and trend data; isolated Ne values are misleading for conservation decisions

## Related Skills

- ecological-genomics/landscape-genomics - Adaptive variation and genotype-environment associations
- ecological-genomics/species-delimitation - Taxonomic unit definition (cite Sukumaran-Knowles caveat for ESU/MU)
- population-genetics/population-structure - Population stratification, STRUCTURE/ADMIXTURE for substructure inference
- population-genetics/selection-statistics - Genome-wide selection signatures including iHS, XP-EHH
- variant-calling/vcf-basics - VCF preparation from RAD-seq or WGS data

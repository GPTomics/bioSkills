# Landscape Genomics Usage Guide

## Overview

Tests genotype-environment associations and identifies loci under local adaptation while correcting for the four-confound landscape (population structure, demographic history, background selection, sampling design). Covers LFMM2 latent factor mixed models with MANDATORY K selection via sNMF cross-entropy elbow (LEA 3), BayPass Core/AUX/C2/IS Bayesian framework with population-covariance matrix Omega (Gautier 2015), RDA / pRDA as the method-of-choice for polygenic adaptation per Forester 2018 (with the genotype-imputation requirement), OutFLANK demographically-calibrated FST outliers, pcadapt PC-based Mahalanobis scans, the Capblancq & Forester 2021 RDA Swiss-army-knife workflow including variance partitioning and adaptive indices, genomic-offset prediction (gradientForest, RONA, LFMM2offset, RDA-offset) WITH the Lind & Lotterhos 2025 three-regime prediction-novelty caveat, Lotterhos & Whitlock 2014 demographic confound + 2015 sampling-design optima (paired contrasts > random > transects under demographic concern), Wang & Bradburd 2014 IBD vs IBE, and Circuitscape + ResistanceGA for landscape resistance.

## Prerequisites

- R with LEA (`BiocManager::install('LEA')`) for sNMF + LFMM2
- pcadapt (`install.packages('pcadapt')`)
- OutFLANK (`devtools::install_github('whitlock/OutFLANK')`)
- vegan for RDA / pRDA (`install.packages('vegan')`)
- gradientForest from R-Forge (`install.packages('gradientForest', repos = 'http://R-Forge.R-project.org')`)
- terra for environmental raster extraction (`install.packages('terra')`)
- qvalue for FDR (`BiocManager::install('qvalue')`)
- BayPass CLI from https://forgemia.inra.fr/mathieu.gautier/baypass_public
- adespatial / Circuitscape (Julia version) / ResistanceGA as needed

## Quick Start

Tell your AI agent what you want to do:

- "Choose K for LFMM2 via sNMF cross-entropy elbow (mandatory; wrong K silently invalidates results)"
- "Run RDA + pRDA per Capblancq & Forester 2021 with imputed genotypes for polygenic adaptation"
- "Use OutFLANK demographically-calibrated FST outliers"
- "Compute BayPass Core (Omega), AUX (env association), C2 (binary contrast)"
- "Predict genomic offset under climate change WITH explicit prediction-novelty regime characterization (Lind 2025)"
- "Optimize sampling design with paired environmental-endpoint contrasts (Lotterhos & Whitlock 2015)"

## Example Prompts

### LFMM2 with K Selection

> "Run sNMF on my VCF for K=1...10 with 5 repetitions per K. Plot cross-entropy and pick K at the elbow. Pass that K to LFMM2 (LEA 3) and report genomic inflation factor (target ~1.0). Re-run at K-1 and K+1 for sensitivity; flag loci detected at only one K as lower-confidence."

### Polygenic RDA: The Capblancq Workflow

> "Implement Capblancq & Forester 2021 RDA workflow: (1) impute missing genotypes BEFORE RDA (mandatory; Forester 2018 benchmark assumes imputed data); (2) forward-select environmental predictors; (3) pRDA conditioning on Q-matrix; (4) variance partitioning into pure-env / pure-spatial / shared / residual; (5) GEA outliers at z > 3 on RDA axes 1-3."

### BayPass Bayesian GEA

> "Convert my VCF to BayPass allele-count format (NOT VCF directly). Run Core to estimate Omega; then AUX for environmental association with Bayes Factor reporting; consider C2 if binary group contrast is the question."

### OutFLANK for Discrete Populations

> "Run OutFLANK with LeftTrim=RightTrim=0.05 and Hmin=0.1 on my 5 discrete populations. Note: OutFLANK requires distinct populations; for continuous gradients use LFMM2 or RDA instead."

### Genomic Offset with the Lind 2025 Caveat

> "Predict genetic offset under SSP5-8.5 2070 using both gradientForest AND RDA-offset. CRITICAL: also compute prediction-novelty per location (Mahalanobis distance to training-data envelope). Annotate offset values with regime: 1 (interpolation), 2 (modest extrapolation), 3 (highly novel; uninformative or misleading per Lind 2025)."

### Landscape Resistance

> "Model landscape resistance with Circuitscape (Julia version) integrating ALL paths. Then optimize the resistance surface with ResistanceGA using MLPE genetic-distance random-effects."

## What the Agent Will Do

1. Identify the four-confound landscape applicable to the dataset (population structure, demographic history, background selection, sampling design) and pick methods that address each
2. For GEA: select K for LFMM2 via sNMF cross-entropy elbow; report sensitivity at K-1, K+1; check genomic inflation factor lambda ~ 1.0
3. For polygenic adaptation: impute missing genotypes BEFORE RDA (mandatory per Forester 2018 benchmark); run pRDA conditioning on Q-matrix; report variance partitioning per Capblancq & Forester 2021
4. For monogenic adaptation in discrete populations: use OutFLANK with trimmed-tail null calibrated by demography
5. For Bayesian inference with population-history correction: BayPass Core + AUX or C2 with Bayes Factor reporting
6. For climate-change prediction: cross-validate genomic offset with gradientForest + RDA-offset + LFMM2offset; characterize prediction-novelty per location; cite Lind & Lotterhos 2025 caveat
7. For sampling design: prefer paired contrasts at environmental endpoints; cite Lotterhos & Whitlock 2015
8. For distinguishing IBD vs IBE: use dbRDA with geographic + environmental matrices (NOT partial Mantel)
9. Apply Storey FDR (q < 0.05) for multiple-testing correction; report number of candidates AND overlap across methods

## Tips

- The number-one failure mode in landscape genomics is treating GEA as a simple signal-extraction problem; ALL FOUR confounds (structure, demography, background selection, sampling) must be addressed
- LFMM2 K selection is the single most consequential parameter; wrong K silently invalidates results (no error message); always use cross-entropy elbow and sensitivity-check at K-1 / K+1
- Genomic inflation factor (GIF / lambda) should be ~1.0; if > 1.5, increase K to absorb more structure; if < 0.5, decrease K
- RDA was established as the polygenic-adaptation gold standard by Forester 2018, but the benchmark assumed IMPUTED genotypes; raw allele matrices with NA degrade RDA performance silently
- Capblancq & Forester 2021 RDA Swiss-army-knife workflow handles variable selection, variance partitioning, GEA outliers, adaptive indices, and genomic offset in one framework
- Storey q-value (`qvalue` package) is more powerful than Benjamini-Hochberg FDR for genomic-scale tests; q < 0.05 is the standard threshold
- OutFLANK requires distinct populations; for continuous gradient sampling, use LFMM2 or RDA instead
- BayPass input is space-separated allele counts, NOT VCF; convert via vcf2baypass.pl or `bcftools query` one-liner
- BayPass Core estimates Omega (the population covariance from genome-wide SNPs); always inspect Omega heatmap before running AUX/C2; it should make biological sense (clustering should reflect known relationships)
- Genomic offset has three regimes per Lind & Lotterhos 2025: interpolation (informative), modest extrapolation (degraded but usable), highly novel (uninformative or misleading); always characterize prediction-novelty
- Sampling design optima per Lotterhos & Whitlock 2015: paired contrasts > random > transects when demographic concern is present (i.e., most landscape-genomics studies)
- LDNe physical-linkage trap: LD-prune SNPs (`PLINK --indep-pairwise 50 5 0.2`) or thin to >= 1 cM apart BEFORE GEA / outlier tests; without LD pruning, outlier counts inflate due to chromosomal clustering
- Gradient forests paper is by Ellis, Smith, Pitcher (2012), NOT by Ellis-Manel-Anderson-Smouse (a common mis-citation in landscape-genomics manuscripts)
- For polyploid species, use polyRAD / fitPoly / updog FIRST to call genotypes correctly; diploid-assuming tools give biased FST and confound GEA
- Circuitscape integrates ALL movement paths weighted by conductance (Julia version is current; Python 4 is deprecated); ResistanceGA optimizes the resistance surface via genetic algorithms
- For distinguishing IBD vs IBE, use dbRDA with geographic and environmental distance matrices; partial Mantel inflates Type I error under spatial autocorrelation (Legendre & Fortin 2010)

## Related Skills

- ecological-genomics/conservation-genetics - Population genetic health (Ne, F_ROH); use ESU/MU framework not species delimitation for management units
- ecological-genomics/community-ecology - Environmental gradient analysis (PERMANOVA + PERMDISP) for species composition
- population-genetics/selection-statistics - Selection scans in human population genetics (iHS, XP-EHH)
- population-genetics/population-structure - STRUCTURE / ADMIXTURE for substructure inference
- variant-calling/vcf-basics - VCF preparation from RAD-seq or WGS

# Landscape Genomics Usage Guide

## Overview

Identifies loci under local adaptation by testing associations between genotypes and environmental variables while correcting for neutral population structure. Covers LFMM2 latent factor mixed models (LEA), pcadapt outlier detection, OutFLANK Fst-based selection scans, RDA-based multi-locus GEA, and gradientForest for predicting genetic vulnerability to climate change. Designed for RADseq, WGS, or SNP array data from organisms sampled across environmental gradients.

## Prerequisites

- R with LEA (`BiocManager::install('LEA')`)
- pcadapt (`install.packages('pcadapt')`)
- OutFLANK (`devtools::install_github('whitlock/OutFLANK')`)
- vegan for RDA-based GEA (`install.packages('vegan')`)
- gradientForest (`install.packages('gradientForest', repos = 'http://R-Forge.R-project.org')`)
- terra for environmental raster extraction (`install.packages('terra')`)
- qvalue for FDR control (`BiocManager::install('qvalue')`)

## Quick Start

Tell your AI agent what you want to do:

- "Find loci associated with temperature in my RADseq data across an elevation gradient"
- "Run pcadapt to detect selection outliers in my SNP dataset"
- "Test genotype-environment associations with LFMM2 controlling for population structure"
- "Predict genetic offset under future climate for my study species"
- "Identify Fst outliers using OutFLANK"

## Example Prompts

### Genotype-Environment Association

> "I have a VCF with 10,000 SNPs from 200 individuals sampled across an elevation gradient. Run LFMM2 to find loci associated with temperature and precipitation, controlling for K=3 population clusters."

> "Test genotype-environment associations using RDA, conditioning on population structure from sNMF. Use WorldClim bioclimatic variables extracted at my sampling coordinates."

### Outlier Detection

> "Run pcadapt on my bed/bim/fam files to detect selection outliers. Help me choose K from the screeplot and filter at q-value < 0.05."

> "Use OutFLANK to identify Fst outlier loci between my 5 sampling populations. Trim the top and bottom 5% of the Fst distribution."

### Climate Vulnerability

> "Predict genetic offset for my study species under SSP5-8.5 2070 climate projections using gradientForest. Show which environmental variables drive the most allele frequency turnover."

### Multi-Method Consensus

> "Run LFMM2, pcadapt, and OutFLANK on my dataset and find the loci identified as candidates by at least two methods."

## What the Agent Will Do

1. Convert genotype data to the required input format (VCF to lfmm/geno/bed)
2. Estimate population structure (K) using sNMF cross-entropy
3. Run the selected GEA method (LFMM2, pcadapt, OutFLANK, RDA) with appropriate parameters
4. Calibrate p-values using genomic inflation factor (LFMM2) or chi-squared fitting (OutFLANK)
5. Apply multiple testing correction (q-value < 0.05 FDR)
6. Report candidate adaptive loci with associated environmental variables
7. Generate diagnostic plots (QQ-plot, Manhattan plot, screeplot, Fst distribution)
8. Optionally predict genetic offset under climate change scenarios with gradientForest

## Tips

- Always estimate K first (sNMF cross-entropy or pcadapt screeplot) before running LFMM2
- Check the genomic inflation factor (GIF/lambda) after LFMM2: target ~1.0; if >1.5, increase K
- pcadapt does not require environmental data, making it useful as a complement to GEA methods
- Combine results from multiple methods (LFMM2 + pcadapt + OutFLANK) for robust candidate lists
- OutFLANK requires distinct populations; continuous sampling designs are better suited to LFMM2 or RDA
- For RDA-based GEA, condition on population structure using the Q-matrix from sNMF as covariates
- Extract environmental data from WorldClim (worldclim.org) using terra::extract() at sampling coordinates
- gradientForest requires at least 10-20 sampling sites for reliable turnover estimation
- Always LD-prune SNPs before pcadapt and OutFLANK to avoid inflated signals from linked loci

## Related Skills

- conservation-genetics - Population genetic health assessment
- community-ecology - Environmental gradient analysis for species
- population-genetics/selection-statistics - Selection scans in human data
- population-genetics/population-structure - Population structure inference
- variant-calling/vcf-basics - VCF input preparation

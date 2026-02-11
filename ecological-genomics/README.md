# ecological-genomics

## Overview

Analyze ecological and environmental genomics data. Covers eDNA metabarcoding, biodiversity metrics, community ecology, genotype-environment associations, conservation genetics, and molecular species delimitation.

**Tool type:** mixed | **Primary tools:** OBITools3, iNEXT, vegan, LEA, hierfstat, ASAP

## Skills
| Skill | Description |
|-------|-------------|
| edna-metabarcoding | Process eDNA amplicon data from raw reads to species occurrence tables |
| biodiversity-metrics | Calculate richness, diversity, and turnover with Hill numbers and iNEXT |
| community-ecology | Analyze community composition with constrained ordination and indicator species |
| landscape-genomics | Test genotype-environment associations and detect local adaptation |
| conservation-genetics | Assess genetic health with Ne estimation, inbreeding, and demographic history |
| species-delimitation | Delimit species from molecular data using distance, tree, and coalescent methods |

## Example Prompts
- "Process my eDNA water samples to identify fish species present"
- "Compare biodiversity across my sampling sites using Hill number rarefaction"
- "Test which environmental variables drive species composition differences"
- "Find loci under local adaptation across an elevation gradient"
- "Estimate effective population size for my endangered species"
- "Delimit species in my COI barcoding dataset"

## Requirements
```bash
# R packages
install.packages(c('vegan', 'iNEXT', 'iNEXT.3D', 'betapart', 'indicspecies',
                    'hierfstat', 'detectRUNS', 'adegenet', 'terra'))
BiocManager::install(c('LEA', 'decontam'))

# R-Forge packages
install.packages('gradientForest', repos = 'http://R-Forge.R-project.org')

# GitHub packages
remotes::install_github('esrud/GONE2')
remotes::install_github('donaldtmcknight/microDecon')
# bPTP (Python package, not R)
# pip install PTP-pyqt5

# Python
pip install OBITools3 cutadapt

# CLI tools
# ASAP: https://bioinfo.mnhn.fr/abi/public/asap/
# Stairway Plot 2: https://github.com/xiaoming-liu/stairway-plot-v2
# NeEstimator: https://github.com/bunop/NeEstimator2.X
```

## Related Skills

- **microbiome** - 16S amplicon processing and clinical microbiome analysis
- **metagenomics** - Shotgun metagenomic classification and profiling
- **population-genetics** - Human population genetics and GWAS
- **phylogenetics** - Tree building for downstream species delimitation
- **comparative-genomics** - Ortholog and synteny analysis across species
- **database-access** - BOLD API and GenBank for reference sequence retrieval
- **variant-calling** - VCF generation from RADseq/WGS for landscape and conservation genomics

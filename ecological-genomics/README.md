# ecological-genomics

## Overview

Analyze ecological and environmental genomics data at postdoc-grade depth. Covers eDNA metabarcoding (ASV vs OTU decision, primer-specific bias, tag-jumping quantification), biodiversity metrics (Hill numbers as effective species counts, coverage-based rarefaction, Chao1 as lower bound), community ecology (mandatory PERMANOVA + PERMDISP, Mantel replacements, JSDMs with residual-covariance-not-interaction caveat), landscape genomics (four-confound GEA framework, RDA for polygenic adaptation with imputation requirement, genomic offset with the Lind 2025 prediction-novelty regime), conservation genetics (Ne by time horizon, the modern 100/1000 rule, Ne/Nc taxonomic variation, genetic-load decomposition), and molecular species delimitation (Sukumaran-Knowles MSC oversplitting critique, DELINEATE remedy, ASAP > ABGD, mPTP > bPTP, integrative taxonomy congruence).

**Tool type:** mixed | **Primary tools:** DADA2, OBITools3, iNEXT, vegan, LEA, BayPass, NeEstimator V2, GONE2, SLiM, ASAP, mPTP, BPP, DELINEATE, Dsuite

## Skills

| Skill | Description |
|-------|-------------|
| edna-metabarcoding | Process eDNA amplicon data with ASV/OTU choice, primer-bias caveats, tag-jumping quantification, decontam screening, and read-counts-not-abundance critique |
| biodiversity-metrics | Compute Hill-number effective-species diversity with coverage-based rarefaction, beta-partition (Baselga + Podani), and SES_MPD with explicit null-model justification |
| community-ecology | Constrained ordination (CCA/RDA/dbRDA) with mandatory PERMANOVA + PERMDISP, JSDMs (HMSC/sjSDM), Mantel replacements, RLQ + fourth-corner |
| landscape-genomics | GEA via LFMM2 (K-elbow), BayPass, RDA polygenic, genomic offset with three-regime caveat, OutFLANK, sampling design optimization |
| conservation-genetics | Ne by time horizon (LDNe / GONE2 / PSMC / dadi / fsc2), F_ROH length-class dating, genetic-load decomposition, 100/1000 rule, tree-sequence SLiM |
| species-delimitation | ASAP (modern successor to ABGD), mPTP, GMYC, BPP-with-calibrated-priors, DELINEATE for Sukumaran-Knowles oversplitting, Dsuite for introgression |

## Example Prompts

- "Process my eDNA water samples to identify fish species with MiFish 12S primers, applying tag-jumping correction (10x heavier for NovaSeq) and decontam as a screening tool"
- "Compare biodiversity across sites using coverage-based Hill-number rarefaction (NOT sample-size rarefaction); report effective species counts, not raw Shannon"
- "Test community-environment relationships with PERMANOVA AND betadisper; report both per Anderson & Walsh 2013"
- "Find polygenic adaptation loci using RDA on imputed genotypes per Capblancq & Forester 2021; predict genomic offset with prediction-novelty regime characterization (Lind 2025)"
- "Estimate Ne by time horizon: GONE2 for recent trajectory (verify cM map column populated); SNeP for RAD-seq physical-linkage correction; PSMC for deep history"
- "Delimit species with ASAP + mPTP, then apply DELINEATE to address MSC oversplitting (Sukumaran-Knowles 2017); validate with integrative taxonomy"

## Requirements

```bash
# R packages
install.packages(c('vegan', 'iNEXT', 'iNEXT.3D', 'betapart', 'indicspecies',
                    'hierfstat', 'detectRUNS', 'adegenet', 'terra', 'picante',
                    'mFD', 'ade4', 'adespatial', 'Hmsc', 'sjSDM', 'qvalue',
                    'pcadapt', 'fossil', 'splits', 'phytools', 'ape',
                    'poppr', 'metabaR', 'occumb'))
BiocManager::install(c('LEA', 'decontam', 'dada2'))

# R-Forge packages
install.packages('gradientForest', repos = 'http://R-Forge.R-project.org')
install.packages('splits', repos = 'http://R-Forge.R-project.org')

# GitHub packages
remotes::install_github('whitlock/OutFLANK')
remotes::install_github('donaldtmcknight/microDecon')
remotes::install_github('esrud/GONE2')

# Python
pip install OBITools3 cutadapt
pip install git+https://github.com/iTaxoTools/PTP-pyqt5  # bPTP

# CLI tools
# ASAP: https://bioinfo.mnhn.fr/abi/public/asap/
# mPTP: https://github.com/Pas-Kapli/mptp (build from C source)
# BPP v4: https://github.com/bpp/bpp
# NeEstimator V2: https://github.com/bunop/NeEstimator2.X (Java)
# SNeP: https://sourceforge.net/projects/snepnetrends/
# GONE2: https://github.com/esrud/GONE2 (CLI form)
# Stairway Plot 2: https://github.com/xiaoming-liu/stairway-plot-v2
# PSMC: https://github.com/lh3/psmc
# SLiM 4: https://messerlab.org/slim/
# Dsuite: https://github.com/millanek/Dsuite
# DELINEATE: https://github.com/jsukumaran/delineate
# Circuitscape: https://github.com/Circuitscape/Circuitscape.jl (Julia)
# BayPass: https://forgemia.inra.fr/mathieu.gautier/baypass_public
```

## Related Skills

- **microbiome** - 16S amplicon processing and clinical microbiome analysis
- **metagenomics** - Shotgun metagenomic classification and profiling
- **population-genetics** - Human population genetics and GWAS
- **phylogenetics** - Tree building for downstream species delimitation
- **comparative-genomics** - Ortholog and synteny analysis across species
- **database-access** - BOLD API and GenBank for reference sequence retrieval
- **variant-calling** - VCF generation from RADseq/WGS for landscape and conservation genomics

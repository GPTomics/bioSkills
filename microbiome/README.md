# microbiome

## Overview

Amplicon (16S rRNA / ITS / 18S) sequencing analysis from demultiplexed reads to compositional differential abundance: DADA2 ASV inference, taxonomy assignment against SILVA/GTDB/Greengenes2/UNITE, alpha/beta diversity, compositional differential abundance, and PICRUSt2 functional prediction. Decision-grade framing throughout: an ASV is a per-run model-inferred sequence not a 97% cluster, 16S resolves to genus not species, microbiome counts are compositional so the DA-tool choice changes the hit list, the diversity answer depends on the rarefaction depth/tree/metric you choose, and PICRUSt2 predicts functional potential not measured activity. It also covers host-organelle (mitochondria/chloroplast) filtering and low-biomass decontamination with decontam. Shotgun metagenomics lives in the metagenomics category.

**Tool type:** mixed | **Primary tools:** DADA2, phyloseq, ALDEx2, QIIME2

## Skills

| Skill | Description |
|-------|-------------|
| amplicon-processing | ASV inference from 16S/ITS amplicon reads with DADA2 (per-run error model, primers-first, truncation budget, chimera removal) |
| taxonomy-assignment | Taxonomic classification of ASVs against SILVA/GTDB/Greengenes2/UNITE (naive Bayes, vsearch-consensus, IDTAXA; region-matched training) |
| diversity-analysis | Alpha/beta diversity with the rarefaction-depth, tree, and metric decisions (phyloseq/vegan + QIIME2 core-metrics; UniFrac, PERMANOVA paired with betadisper) |
| differential-abundance | Compositional differential abundance via a consensus of ALDEx2, ANCOM-BC2, MaAsLin2, and LinDA (Nearing-2022 tool disagreement) |
| functional-prediction | Predict functional potential from 16S with PICRUSt2 (NSTI-gated, potential not activity, taxonomy-conditioned) |
| qiime2-workflow | QIIME2 artifact/provenance framework that ties the amplicon pipeline together (.qza/.qzv, semantic types, import/export, provenance replay) |

## Example Prompts

- "Process paired-end 16S reads and infer ASVs with a per-run DADA2 error model"
- "Assign taxonomy to my ASVs against SILVA 138 with a region-matched classifier"
- "Calculate alpha and beta diversity, choose a sampling depth, and compare groups with PERMANOVA"
- "Find differentially abundant taxa with a consensus of compositional methods"
- "Predict KEGG/MetaCyc functional potential from 16S with PICRUSt2 and report NSTI"
- "Filter host mitochondria/chloroplast reads and remove reagent contaminants from low-biomass samples"

## Requirements

```bash
# CLI (bioconda)
conda install -c bioconda cutadapt itsxpress

# QIIME2 installs as its own conda env; verify the current amplicon distribution
# installer at https://docs.qiime2.org

# PICRUSt2 (own conda env)
conda install -c bioconda picrust2
```

```r
# R / Bioconductor
BiocManager::install(c('dada2', 'DECIPHER', 'phyloseq', 'ALDEx2', 'ANCOMBC', 'Maaslin2'))
install.packages(c('vegan', 'picante'))
```

## Related Skills

- **metagenomics** - Shotgun (whole-genome) taxonomic and functional profiling; shared compositional and diversity theory
- **pathway-analysis** - Enrichment analysis of predicted KEGG/MetaCyc functions
- **data-visualization** - Ordination plots, taxonomic barplots, heatmaps
- **read-qc** - Adapter/primer trimming and read quality reports before DADA2
- **phylogenetics** - Phylogenetic tree I/O for UniFrac and Faith PD

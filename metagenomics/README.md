# metagenomics

## Overview

Read-based taxonomic and functional profiling of shotgun metagenomes. Decision-grade framing: a metagenomic result is a position in a choice-chain (extraction -> depletion -> depth -> classifier -> database -> normalization), never a direct observation - so the database defines what can be detected, read counts are not abundances, the data are compositional, and absence usually means not-detectable-by-this-chain. Assembly and MAG recovery live in genome-assembly/metagenome-assembly.

**Tool type:** cli | **Primary tools:** Kraken2, MetaPhlAn, Bracken, HUMAnN, AMRFinderPlus, inStrain, decontam

## Skills

| Skill | Description |
|-------|-------------|
| contamination-controls | Host depletion, kitome blanks/decontam, mocks, and depth checks before profiling |
| kraken-classification | K-mer classification with Kraken2 + Bracken; confidence and false-positive control |
| metaphlan-profiling | Marker-gene species/SGB profiling; cell fraction, not read fraction |
| abundance-estimation | Bracken re-estimation plus compositional, normalization, and absolute-load handling |
| functional-profiling | HUMAnN gene-family and pathway abundance; potential, not activity |
| amr-detection | Community resistome from reads/contigs; an ARG hit is not a phenotype |
| strain-tracking | inStrain/StrainPhlAn strain resolution; a strain is a threshold, not a thing |
| metagenome-visualization | Honest figures and community statistics from profiler tables |

## Example Prompts

- "Remove human reads and flag reagent contaminants before profiling my low-biomass samples"
- "Classify my shotgun reads at a confidence threshold, then estimate species-level abundance"
- "Profile species with a marker-gene method and keep the unknown fraction for my environmental samples"
- "CLR-transform my abundance table and run a compositional differential-abundance test"
- "Quantify the resistome of my metagenome reads and report ARG abundance, not resistance"
- "Detect whether a strain is shared between my paired samples"
- "Profile functional pathways without dropping the unmapped fraction"

## Requirements

```bash
# Profiling
conda install -c bioconda kraken2 bracken metaphlan humann

# Resistome, strains, controls
conda install -c bioconda ncbi-amrfinderplus rgi instrain drep skani hostile nonpareil krona
Rscript -e "BiocManager::install(c('decontam','phyloseq','ALDEx2'))"
```

## Related Skills

- **genome-assembly** - Assemble the community into MAGs (the assembly-based view; this category is read-based)
- **read-qc** - Adapter/quality trimming and host screening before profiling
- **epidemiological-genomics** - Isolate-level AMR surveillance and outbreak typing from pure cultures
- **pathway-analysis** - Organism-centric functional enrichment downstream of HUMAnN
- **microbiome** - Amplicon/16S/QIIME2 diversity and differential abundance

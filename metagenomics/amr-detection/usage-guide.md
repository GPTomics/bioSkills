# AMR Detection (Community Resistome) - Usage Guide

## Overview
This skill profiles the antimicrobial-resistance gene content of shotgun metagenomes: read-based quantification (RGI bwt, AMR++/MEGARes, ARGs-OAP/SARG, deepARG, GROOT) and presence calling on assembled contigs or MAGs (AMRFinderPlus, ABRicate). An ARG hit is a sequence match, not a phenotype - in a mixed community it has no host and no genomic context until you assemble, and assembly tends to break exactly where the ARGs are. The honest output is "ARG present at abundance X," never "the sample is resistant." Pure-culture isolate AMR, point mutations, and phenotype/MIC prediction belong to epidemiological-genomics/amr-surveillance.

## Prerequisites
```bash
conda install -c bioconda ncbi-amrfinderplus rgi abricate
amrfinder -u            # update the NCBI Reference Gene Catalog
rgi load --card_json /path/to/card.json --local   # load CARD for RGI
abricate --setupdb      # or abricate-get_db
```

Conceptual prerequisites:
- Decide read-based (quantitative, no host) vs assembly-based (presence per contig, possible host) up front.
- Gene fraction / breadth-of-coverage is mandatory for read-based calls; identity alone counts fragments as genes.
- Point mutations need a single known organism - largely out of reach for a mixed community.
- Cross-study ARG abundance is rarely comparable; record DB version, tool version, normalization unit, and depth.

## Quick Start
Tell your AI agent what you want to do:
- "Quantify the resistome of my metagenome reads with RGI bwt and report ARG abundance, not resistance"
- "Call ARGs on my assembled MAGs with AMRFinderPlus"
- "Normalize ARG abundance to copies per 16S across my samples"
- "Figure out whether this carbapenemase is on a plasmid or in a pathogen"

## Example Prompts

### Read-based resistome
> "I have 30 sewage shotgun samples. Quantify the resistome from reads with a gene-fraction filter, normalize per 16S, and compare relative ARG abundance between sites - within this one pipeline only."

### Contig/MAG presence
> "Call ARGs on my binned MAGs with AMRFinderPlus, and tell me which MAG carries each gene. Do not call point mutations unless the MAG is a single resolved species."

### Host/MGE linkage
> "A blaNDM read mapped in my metagenome. Help me determine which organism carries it and whether it is mobile."

### Avoiding the phenotype trap
> "My pipeline output says 'resistant to colistin.' Reframe this correctly for a metagenome."

## What the Agent Will Do
1. Choose read-based or assembly-based based on whether quantification, host context, or point mutations are needed.
2. Apply a gene-fraction/breadth filter to read-based calls and a curated per-gene cutoff to contig calls.
3. Report ARG presence and relative abundance with the normalization unit, never a resistance phenotype.
4. Route host/MGE-linkage questions to assembly+binning, long reads, or Hi-C.
5. Defer pure-culture isolate AMR, point mutations, and MIC prediction to epidemiological-genomics/amr-surveillance.

## Tips
- AMRFinderPlus uses curated per-gene cutoffs; do not override `--ident_min` with a global value.
- ABRicate is contigs-only and acquired-genes-only - it never reports a point mutation.
- RGI bwt cannot screen the resistance SNP on variant models; a gyrA read hit is not a resistance call.
- Record `amrfinder -V` and `abricate --list` versions; database snapshots change calls.
- Compare ARG abundance only within a study, same DB/normalization/depth.

## Read vs Assembly

| Axis | Read-based | Assembly-based |
|------|-----------|----------------|
| Low-abundance sensitivity | high | low |
| Quantification | yes | presence/absence |
| Host / MGE context | none without binning | possible |
| Point mutations | unreliable | yes, with organism/model |

## Resources
- [AMRFinderPlus wiki](https://github.com/ncbi/amr/wiki)
- [CARD / RGI](https://card.mcmaster.ca/)
- [MEGARes / AMR++](https://www.meglab.org/)

## Related Skills

- epidemiological-genomics/amr-surveillance - Isolate AMR, point mutations, phenotype/MIC, typing
- functional-profiling - General gene-family/pathway abundance
- kraken-classification - Taxonomic context for the community
- genome-assembly/metagenome-assembly - Assembly/binning and ARG-host linkage
- contamination-controls - Host depletion before resistome profiling
- workflows/metagenomics-pipeline - End-to-end shotgun analysis

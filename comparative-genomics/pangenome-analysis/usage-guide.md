# Pangenome Analysis - Usage Guide

## Overview

The pangenome is the union of all genes across a sampled group of genomes. The Tettelin partition (Tettelin 2005 PNAS 102:13950) splits it into core (universal), shell (in many), cloud (rare), and species-specific genes. **Bacterial pangenome** tools (Panaroo, PPanGGOLiN, PEPPAN) cluster genes into orthogroups for compact genomes. **Eukaryotic pangenome** tools (Minigraph-Cactus, PGGB, vg) build graph-based representations of haplotype-level variation. The Roary tool (Page 2015) is deprecated; Panaroo (Tonkin-Hill 2020) reduces accessory inflation by 30-50% on the same data.

For **repetitive and clinically-relevant genes** (MHC class II, DAZ1-4, OPN1LW/OPN1MW), PGR-TK (Chin 2023 Nat Methods 20:1213) uses MAP graphs + principal bundle decomposition; its successor PANGEA is in development.

## Prerequisites

```bash
# Bacterial pangenome
conda install -c bioconda panaroo ppanggolin peppan get_homologues anvio

# Eukaryotic pangenome
conda install -c bioconda cactus pggb vg pangenie

# PGR-TK / PANGEA
conda install -c bioconda pgr-tk
# Or Docker: quay.io/cschin/pgr-tk

# Annotation
conda install -c bioconda bakta prokka

# Pan-GWAS
conda install -c bioconda scoary pyseer

# Recombination
conda install -c bioconda clonalframeml

# QC
conda install -c bioconda busco compleasm checkm2
```

## Quick Start

Tell the AI agent what to do:
- "Build a bacterial pangenome from these 100 E. coli genomes using Panaroo, then partition into core/shell/cloud with PPanGGOLiN"
- "Construct a human pangenome graph from 50 haplotype-resolved assemblies using Minigraph-Cactus with GRCh38 reference"
- "Apply PGR-TK to analyze the MHC class II locus across diverse human pangenome haplotypes, producing principal bundle decomposition for the region"
- "Run pan-GWAS with Scoary to identify genes associated with antibiotic resistance phenotype across the pangenome"

## Example Prompts

### Bacterial Pangenome

> "Build a high-quality pangenome from 200 Mycobacterium tuberculosis genomes. First re-annotate all genomes with Bakta. Then run Panaroo in strict mode for clonal pangenomes. Apply Heaps law fitting on resampling to determine if the pangenome is open or closed. Report core, shell, and cloud gene counts."

### Eukaryotic Pangenome Graph

> "Construct a pangenome graph from 50 plant haplotype-resolved assemblies (Arabidopsis thaliana strains and Arabidopsis lyrata) using Minigraph-Cactus with the TAIR10 reference. Produce GFA, VCF, and GBZ outputs. Genotype short-read data from a held-out sample against this graph using vg giraffe and PanGenie."

### Clinical Gene Analysis with PGR-TK

> "Apply PGR-TK to analyze the HLA-DRB1 locus across the HPRC pangenome (90 haplotypes). Generate principal bundle decomposition for the region; report the dominant haplotype bundles. Note: PGR-TK is being succeeded by PANGEA (developed by DGI / Diploid Genomics); check the cschin/pgr-tk upstream repository for the latest pointer."

## What the Agent Will Do

1. **Validate inputs**: BUSCO/Compleasm completeness per genome; CheckM2 for MAGs; consistent annotation pipeline
2. **Pre-annotate consistently** with Bakta (bacteria) or projection from reference (eukaryotes)
3. **Run pangenome analysis** with appropriate tool:
   - Bacterial: Panaroo (strict mode) + PPanGGOLiN (partition)
   - Eukaryotic: Minigraph-Cactus (production) or PGGB (reference-free)
   - Repetitive / clinical: PGR-TK with principal bundle decomposition
4. **Apply Tettelin partition** at 95% / 99% core thresholds
5. **Fit Heaps law** for open / closed pangenome characterization
6. **Cross-validate** with PPanGGOLiN HMM partition for principled class boundaries
7. **Pan-GWAS** with Scoary or pyseer if phenotype available
8. **Report**: core/shell/cloud sizes, openness alpha, per-genome accessory complement
9. **Caveats**: annotation heterogeneity, assembly fragmentation, recombination

## Tips

- Use Bakta (not Prokka) for current bacterial annotation; GenBank-compliant + faster
- Panaroo `--clean-mode strict --remove-invalid-genes` reduces annotation-error inflation 30-50%
- Roary is deprecated; migrate to Panaroo
- For clonal lineages (Mtb), Panaroo strict mode is gold standard
- For genus-scale (1000+ genomes), PEPPAN scales better than Panaroo
- PPanGGOLiN HMM partition (persistent/shell/cloud) more principled than Tettelin threshold
- Heaps law requires >= 100 genomes; report 95% CI from resampling
- Eukaryotic pangenome: Minigraph-Cactus for production; PGGB for reference-free
- PGR-TK is purpose-built for repetitive / clinically-relevant genes (MHC, DAZ1-4, OPN1LW); use it when standard tools collapse alleles
- PGR-TK is archived (April 2026); PANGEA is the successor, developed by DGI (Diploid Genomics); check upstream repo for current pointer
- For variant calling from pangenome graph, use vg giraffe + vg call
- For SV genotyping, use PanGenie
- Pre-filter MAGs with CheckM2: >= 70% completeness, < 5% contamination
- Recombination must be addressed for highly recombinogenic species (Neisseria, S. pneumoniae)
- Document annotation pipeline + version; mixing pipelines inflates accessory genome

## Related Skills

comparative-genomics/whole-genome-alignment - Minigraph-Cactus / PGGB use Cactus
comparative-genomics/ortholog-inference - Bacterial pangenome clusters are orthogroups
comparative-genomics/hgt-detection - Accessory genes often HGT-derived
comparative-genomics/gene-family-evolution - CAFE5 on Panaroo presence/absence
comparative-genomics/genome-distance-and-species-delineation - ANI for clustering precedes pangenome
genome-annotation/prokaryotic-annotation - Bakta input
genome-annotation/repeat-annotation - Repeat masking for eukaryotic pangenome
variant-calling/structural-variant-calling - Pangenome SV calling
metagenomics/amr-detection - AMR genes in bacterial accessory

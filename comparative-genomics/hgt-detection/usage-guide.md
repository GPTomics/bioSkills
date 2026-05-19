# HGT Detection - Usage Guide

## Overview

Horizontal gene transfer (HGT / LGT) detection combines three orthogonal signal classes: **composition** (GC%, codon usage, tetranucleotide z-scores -- recent transfers), **phylogeny** (gene tree-species tree incongruence -- any age if signal preserved), and **comparative** (orthogroup presence/absence -- patchy distribution). The 2024 standard is **ALE / GeneRax / AleRax** for probabilistic DTL reconciliation, complemented by **HGTector v2** (BLAST distribution) and **AvP** / **HGTphyloDetect** for eukaryotic analysis. **Contamination is the dominant false-positive source in eukaryotic HGT analysis**: the Boothby 2015 tardigrade "17% HGT" claim was debunked (Koutsovoulos 2016) and Crisp 2015 "145 HGTs in humans" was reduced ~90% after Salzberg 2017 contamination-aware reanalysis. FCS-GX or BlobTools contamination filter is mandatory before eukaryotic HGT calls.

Differential gene loss (DGL) mimics HGT signatures; explicit testing via ALE event posteriors required.

## Prerequisites

```bash
# Probabilistic DTL
conda install -c bioconda ale generax
git clone --recursive https://github.com/BenoitMorel/AleRax && cd AleRax && ./install.sh

# BLAST-distribution screening
pip install hgtector

# Eukaryotic HGT pipeline (AvP)
git clone https://github.com/GDKO/AvP && cd AvP && pip install .

# Mobile elements (bacterial)
# IslandViewer 4 (web): https://www.pathogenomics.sfu.ca/islandviewer/
git clone https://github.com/clb21565/mobileOG-db

# Contamination filtering (MANDATORY for eukaryotic HGT)
conda install -c bioconda fcs-gx blobtoolkit kraken2

# AMR + mobile element annotation
conda install -c bioconda amrfinder

# Tree inference for phylogenetic HGT
conda install -c bioconda iqtree
```

## Quick Start

Tell the AI agent what to do:
- "Screen this bacterial genome for HGT using HGTector v2 + IslandViewer 4"
- "Run ALE on UFBoot gene trees against species tree to identify per-branch transfer events"
- "Apply AvP for phylogenetic HGT detection in a eukaryotic genome (with prior contamination filter)"
- "Compute amelioration timeline for candidate HGTs to estimate transfer age"

## Example Prompts

### Bacterial Genome HGT Screen

> "Screen 50 bacterial genomes for HGT using HGTector v2 + ALE undated reconciliation. Filter contamination with FCS-GX first. Annotate candidate transfers with mobileOG-db (mobile elements) and AMRFinderPlus (resistance). Report per-genome HGT count and per-branch transfer rate."

### Eukaryotic HGT (Contamination-Cleared)

> "Identify HGT in a recently sequenced bdelloid rotifer genome. MANDATORY first step: FCS-GX contamination filter, then BlobTools coverage-GC visualization. Then run AvP for phylogenetic HGT detection. Cross-validate with composition (GC z-score, codon usage). Report only HGT candidates that pass contamination filter + show concordance across composition + phylogeny."

### Distinguish HGT from Differential Gene Loss

> "Two related bacterial genera share a gene that's absent in their immediate sister lineage. Test via ALE undated whether this is HGT or DGL. Inspect branch-wise event posteriors (T vs L); if loss rate at sister branches is non-negligible, prefer the loss explanation. Cross-check with RANGER-DTL sensitivity sweep across cost weights."

## What the Agent Will Do

1. **Validate inputs**: assembly quality (BUSCO, N50); contamination check for eukaryotes (FCS-GX, BlobTools, Kraken2)
2. **For bacterial single genome**: IslandViewer 4 + HGTector v2
3. **For 5-200 bacterial genomes**: ALE undated on UFBoot orthogroup gene trees
4. **For eukaryotic**: AvP / HGTphyloDetect after contamination cleanup
5. **Compute composition**: GC% z-scores; codon adaptation index; tetranucleotide bias
6. **Run phylogenetic incongruence**: AU/SH topology tests for individual gene trees
7. **Apply concordance filter**: >= 2 of 3 signal classes
8. **Annotate mobile elements**: integrases, transposases, phage proteins
9. **Distinguish HGT from DGL**: ALE event posteriors (T vs L on candidate branch)
10. **Report**: candidate HGTs with multi-class evidence, donor lineage, transfer age
11. **Caveats**: contamination, amelioration, ILS confound, gBGC

## Tips

- ALE / GeneRax / AleRax have largely superseded parsimony reconciliation for serious work
- Contamination filter is MANDATORY before any eukaryotic HGT analysis (Boothby / Crisp lessons)
- FCS-GX is now required by NCBI GenBank for new eukaryote submissions
- Amelioration (Lawrence-Ochman 1998) erodes compositional signal over 50-200 Myr
- Use phylogenetic methods for HGT older than ~50 Myr
- Concordance across 2+ signal classes is required; single-method claims are weak
- HGTector v2 for fast bacterial screen; ALE for posterior support on candidates
- AvP for eukaryotic; orthogroup-aware end-to-end pipeline
- DGL mimics HGT; inspect ALE T vs L posteriors at sister branches
- gBGC confounds compositional analyses; check W->S substitution ratios
- For close-relative donors, BLAST-distribution methods (HGTector / DarkHorse / AI) fail; use phylogenetic methods
- Donor lineage requires >= 5 sampled genomes from candidate order
- Mobile-element signatures (integrase, transposase, phage, tRNA) corroborate transfer
- Apply AMRFinderPlus for resistance genes on putative HGT islands

## Related Skills

comparative-genomics/gene-tree-species-tree-reconciliation - Deep DTL methodology
comparative-genomics/ortholog-inference - Orthogroup construction for ALE input
comparative-genomics/introgression-detection - Distinguish HGT from hybridization
comparative-genomics/genome-distance-and-species-delineation - ANI / dDDH for donor distance
comparative-genomics/pangenome-analysis - Pangenome accessory genes
phylogenetics/modern-tree-inference - Gene-tree inference with UFBoot
phylogenetics/bayesian-inference - PhyloBayes-MPI CAT-GTR for deep-node LBA
alignment/multiple-alignment - Codon-aware alignment
metagenomics/amr-detection - Mobile element / resistance context
genome-annotation/prokaryotic-annotation - Annotation alongside HGT detection

# Synteny Analysis - Usage Guide

## Overview

Synteny analysis detects conserved gene order across genomes and identifies structural rearrangements. **MCScanX** (Wang 2012) is the most widely used dynamic-programming tool for collinear gene blocks. **JCVI / MCScan Python** (Tang 2008) provides publication-grade visualization. **GENESPACE** (Lovell 2022) is the modern standard for plant comparative genomics with riparian plots and pan-gene tracks. **SyRI** (Goel 2019) handles structural variant detection from whole-genome alignment.

The fundamental distinction: synteny != collinearity (synteny is "same chromosome"; collinearity is "same order"). The choice of tool depends on whether the question is gene-level co-linearity (MCScanX, JCVI), whole-genome structural rearrangements (SyRI, AnchorWave), multi-genome macrosynteny (GENESPACE, ntSynt), or synteny-aware orthology (GENESPACE, ProteinOrtho-synteny).

For whole-genome alignment substrate (Cactus / Minigraph-Cactus), see [[whole-genome-alignment]]. For WGD-aware analysis (Ks plots), see [[whole-genome-duplication]]. For pangenome graphs, see [[pangenome-analysis]].

## Prerequisites

```bash
# MCScanX (compile from source)
git clone https://github.com/wyp1125/MCScanX && cd MCScanX && make

# JCVI / MCScan Python
pip install jcvi

# GENESPACE (R)
Rscript -e "remotes::install_github('jtlovell/GENESPACE')"
# GENESPACE bundles OrthoFinder 2.5.x specifically

# SyRI + plotsr
conda install -c bioconda syri plotsr

# AnchorWave (WGD-aware)
conda install -c bioconda anchorwave

# i-ADHoRe 3.0 (deeply diverged species)
git clone https://github.com/VIB-PSB/i-ADHoRe && cd i-ADHoRe && mkdir build && cd build && cmake .. && make

# ntSynt (multi-genome macrosynteny)
pip install ntsynt

# Pairwise WGA tools
conda install -c bioconda minimap2 mummer4

# Repeat masking (MANDATORY before synteny)
conda install -c bioconda repeatmodeler repeatmasker
```

## Quick Start

Tell the AI agent what to do:
- "Detect syntenic blocks between two plant genomes using MCScanX"
- "Build a riparian plot of synteny across 10 grass genomes with GENESPACE"
- "Find structural rearrangements (inversions, translocations, duplications) with SyRI"
- "Generate publication-grade dotplot with JCVI for human vs mouse synteny"

## Example Prompts

### Plant Comparative Genomics

> "Compare Arabidopsis thaliana, Brassica napus, Camelina sativa, and Capsella rubella using GENESPACE. Build the riparian plot with orthogroup-constrained anchors. Identify lineage-specific gene gains and losses via pan-gene tracks. Cross-reference with [[whole-genome-duplication]] for WGD context."

### Structural Variant Inventory

> "Identify INVs, TRANSs, and DUPs between two Solanum lycopersicum cultivars using minimap2 -x asm5 + SyRI + plotsr. Filter SyRI INVs to > 5 kb (biological-relevance threshold). Cross-validate with PacBio long-read SV calls."

### Microsynteny Network Analysis

> "Build a synteny network across 50 angiosperm genomes using SynNet. Identify microsynteny-conserved gene clusters (metabolic operons, NLR clusters) that persist across deep divergence. Compare with macrosynteny decay (Murat 2010)."

## What the Agent Will Do

1. **Validate inputs**: assembly N50 (>= 1 Mb); BUSCO/Compleasm completeness; repeat masking
2. **Pre-mask repeats**: RepeatModeler2 species-specific library + RepeatMasker softmask
3. **All-vs-all BLAST**: DIAMOND (faster) or BLASTP+
4. **Run synteny detection** with appropriate tool:
   - MCScanX for general use
   - GENESPACE for plant comparative
   - JCVI for publication plots
   - SyRI for structural variants
   - AnchorWave for WGD lineages
5. **Tandem duplicate collapse**: automatic in MCScanX (5-gene window default)
6. **Classify synteny relationships**: 1:1, 1:many, many:many
7. **For SVs**: minimap2 + SyRI; filter INVs > 5 kb
8. **Visualize**: dotplots, karyotypes, riparian plots, plotsr SV plots
9. **Cross-validate**: MCScanX vs JCVI (consistent stringency)
10. **Report**: block-level statistics, synteny scale (micro vs macro), per-chromosome relationships
11. **Caveats**: repeat-driven false synteny, fragmentation, tandem inflation, reference-guided bias, WGD effects

## Tips

- Soft-mask repeats before BLAST; unmasked TEs produce millions of paralogous hits
- Minimum N50 of 1 Mb for reliable synteny detection
- Reference-guided scaffolding creates circular reasoning if used for self-comparison
- MCScanX block size minimum: -s 5 (5 anchor genes); raise to 10 for stringent
- GENESPACE is the modern standard for plant comparative; OrthoFinder + MCScanX-constrained
- For animal/fungal multi-genome, JCVI is most generally useful
- SyRI requires chromosome-level assemblies; works on pairs
- AnchorWave WGD-aware proali mode requires explicit ploidy
- For deeply diverged species (< 70% identity), i-ADHoRe or LASTZ chains/nets
- minimap2 -x asm5 for closely related; LASTZ for distant
- ntSynt for alignment-free multi-genome macrosynteny
- For Ks-dating of synteny blocks, see [[whole-genome-duplication]]
- Tandem duplicates collapsed by MCScanX 5-gene window default; tune for NLR clusters
- Microsynteny vs macrosynteny: minimum 5 anchors microsynteny; 20+ macrosynteny
- For WGA substrate at clade level, see [[whole-genome-alignment]]
- Document microsynteny vs macrosynteny scale in any synteny claim

## Related Skills

comparative-genomics/whole-genome-duplication - Ks-based WGD detection, KsRates
comparative-genomics/whole-genome-alignment - Cactus / Minigraph-Cactus / LASTZ
comparative-genomics/ortholog-inference - GENESPACE depends on OrthoFinder
comparative-genomics/pangenome-analysis - Bacterial pangenome (different problem)
comparative-genomics/comparative-annotation-projection - TOGA / CESAR via WGA-derived synteny
comparative-genomics/positive-selection - Selection on syntenic gene pairs
phylogenetics/modern-tree-inference - Phylogenetic context for dating
genome-assembly/assembly-qc - BUSCO before synteny
alignment/structural-alignment - Sequence-level underlying SyRI
variant-calling/structural-variant-calling - Short-read SV (orthogonal to SyRI WGA)

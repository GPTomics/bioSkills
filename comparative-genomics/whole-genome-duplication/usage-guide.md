# Whole Genome Duplication - Usage Guide

## Overview

Whole-genome duplication (WGD / paleopolyploidy) events shape genome evolution by doubling all genes; subsequent gene loss is biased toward certain functional categories (Birchler & Veitia 2007 Plant Cell 19:395 gene balance hypothesis). Detection combines **Ks distribution** (synonymous substitution rates among paralog pairs) with **synteny block analysis** (parallel collinear blocks). The 2024 standard pipeline is **wgd v2** (Chen & Zwaenepoel 2024 Bioinformatics 40:btae272). For comparing WGDs across lineages, **KsRates** (Sensalari 2022 Bioinformatics 38:530) is mandatory because substitution rates vary across the tree.

The vertebrate 2R, teleost 3R, and salmonid Ss4R WGDs are well-established; plant WGDs are extremely common (every angiosperm has at least one ancestral WGD).

## Prerequisites

```bash
# wgd v2 + KsRates
pip install wgd ksrates

# DupGen_finder (Perl)
git clone https://github.com/qiao-xin/DupGen_finder

# MCScanX (dependency for wgd syn)
git clone https://github.com/wyp1125/MCScanX && cd MCScanX && make

# Phylogenomic WGD detection
git clone https://bitbucket.org/barker-lab/maps  # MAPS
julia -e 'using Pkg; Pkg.add("Whale")'           # Whale.jl

# PAML for Ks computation
conda install -c bioconda paml

# R packages for mixture models
Rscript -e "install.packages(c('mclust', 'rmixmod'))"
```

## Quick Start

Tell the AI agent what to do:
- "Detect whole-genome duplications in this plant genome using wgd v2 Ks distribution analysis"
- "Position the WGD event relative to speciation using KsRates with three outgroup species"
- "Classify all paralog pairs as tandem, proximal, segmental, or WGD using DupGen_finder"
- "Apply Whale.jl Bayesian DTL+WGD modeling with explicit WGD nodes on the species tree"

## Example Prompts

### Plant WGD Discovery

> "Detect WGD events in the Arabidopsis thaliana genome (CDS + protein + GFF provided). Build the paranome with wgd dmd, compute the Ks distribution with wgd ksd, identify synteny anchors with wgd syn, and fit Gaussian mixture models to determine the number of WGD events. Report the inferred age of each WGD peak."

### Cross-Lineage WGD Comparison

> "Compare WGD events in Brassica napus to Arabidopsis thaliana, correcting for substitution rate variation. Use KsRates with at least two outgroup species (Carica papaya and Vitis vinifera). Position the WGD events relative to known speciation events on the Brassicales tree. Determine whether B. napus's WGD predates or postdates the divergence from A. thaliana."

### Phylogenomic WGD Placement

> "Place the WGD events on the vertebrate species tree using Whale.jl with native WGD modeling. The species tree should have explicit WGD nodes at the gnathostome (2R) and teleost (3R) ancestors. Report posterior probabilities for each WGD and retention parameters for each. Cross-validate with MAPS phylogenomic placement."

## What the Agent Will Do

1. **Validate inputs**: Compleasm completeness; CDS sequence validity; species tree consistency
2. **Build paranome** with wgd dmd (DIAMOND-based)
3. **Compute Ks distribution** with wgd ksd (codon-aware alignment + PAML yn00 / codeml)
4. **Identify synteny anchors** with wgd syn (MCScanX-based)
5. **Apply KsRates** for substitution rate correction (if cross-lineage comparison)
6. **Classify paralog pairs** with DupGen_finder (tandem / proximal / dispersed / segmental / WGD)
7. **Fit mixture models** (GMM, ELMM) with BIC-based model selection
8. **Report**: WGD count, Ks peak locations + ages, synteny block synchrony, KsRates-corrected positioning
9. **Cross-validate** with MAPS / Whale.jl for phylogenomic placement
10. **Caveats**: Ks saturation; tandem-driven peaks; substitution rate variation

## Tips

- wgd v2 (heche-psb/wgd) replaces deprecated wgd v1 (arzwa/wgd) -- migrate any old scripts
- KsRates is mandatory for cross-lineage WGD comparison; rates vary 2-5x across angiosperms
- Restrict Ks-based analysis to Ks < 1.5; > 2 saturated and unreliable
- Use DupGen_finder to separate tandem-driven peaks from true WGD peaks
- Synteny block synchrony (per-block Ks IQR overlap) confirms true WGD over sequential duplications
- BIC-based mixture-model selection; report 1-5 components
- For polyploids, assign subgenomes before WGD analysis (SubPhaser, GENESPACE)
- For vertebrate 2R and teleost 3R, Ks is saturated; use MAPS / Whale.jl phylogenomic placement
- Salmonid Ss4R (80-100 Myr) is recent enough for Ks-based dating
- For new WGD claims, require concordance across Ks + synteny + KsRates + phylogenomic placement
- BUSCO completeness > 90% before running wgd v2
- Document outgroup choice for KsRates (>= 2 outgroups at different distances)
- Whale.jl native WGD modeling preferred over interpreting D burst from ALE

## Related Skills

comparative-genomics/synteny-analysis - Synteny-anchored Ks (wgd syn / GENESPACE)
comparative-genomics/ortholog-inference - Ortholog detection feeds DupGen / wgd
comparative-genomics/gene-tree-species-tree-reconciliation - Whale.jl native WGD modeling
comparative-genomics/gene-family-evolution - CAFE5 birth-death across families
comparative-genomics/positive-selection - Selection on retained ohnologs post-WGD
phylogenetics/divergence-dating - Time-calibrated WGD positioning
alignment/multiple-alignment - Codon-aware MSA for Ks
genome-assembly/assembly-qc - BUSCO before WGD analysis

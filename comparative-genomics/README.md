# comparative-genomics

## Overview

Comparative genomics for analyzing genome evolution across species: synteny and collinearity, positive selection detection, ancestral sequence reconstruction, ortholog inference, and horizontal gene transfer detection.

**Tool type:** mixed | **Primary tools:** MCScanX, PAML, OrthoFinder, HyPhy, SyRI

## Skills

| Skill | Description |
|-------|-------------|
| synteny-analysis | Genome collinearity, WGD detection, structural rearrangements with MCScanX, SyRI, JCVI |
| positive-selection | dN/dS site/branch/branch-site tests with PAML codeml and HyPhy (BUSTED, MEME, FEL, FUBAR) |
| ancestral-reconstruction | Ancestral sequence reconstruction with PAML and IQ-TREE for protein resurrection |
| ortholog-inference | Ortholog/paralog classification with OrthoFinder (tree-based) and ProteinOrtho (graph-based) |
| hgt-detection | Horizontal gene transfer via compositional, phylogenetic, and taxonomic distribution methods |

## Example Prompts

- "Find syntenic blocks between human and mouse genomes and check for WGD signatures"
- "Test for positive selection on this gene across mammals using PAML and cross-validate with HyPhy MEME"
- "Reconstruct the ancestral sequence at the mammalian root for protein resurrection"
- "Find orthologs of BRCA1 across vertebrates and classify by copy number"
- "Detect horizontally transferred genes in this bacterial genome using compositional and phylogenetic methods"
- "Screen this gene family for recombination with GARD, then run branch-site selection tests"

## Requirements

```bash
# OrthoFinder (tree-based orthology)
conda install -c bioconda orthofinder

# PAML (selection analysis, ASR)
conda install -c bioconda paml

# HyPhy (modern selection testing: BUSTED, MEME, FEL, GARD)
conda install -c bioconda hyphy

# PRANK (recommended aligner for selection studies)
conda install -c bioconda prank

# MCScanX (synteny detection)
git clone https://github.com/wyp1125/MCScanX
cd MCScanX && make

# SyRI (structural rearrangements)
pip install syri

# JCVI (synteny visualization)
pip install jcvi
```

## Related Skills

- **phylogenetics** - Tree building for selection analysis and ASR
- **alignment** - MSA preparation for codon analysis (PRANK recommended)
- **population-genetics** - Within-species selection statistics
- **genome-annotation** - Annotation transfer via orthology and synteny

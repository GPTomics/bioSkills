# Ortholog Inference - Usage Guide

## Overview

Ortholog inference distinguishes genes whose history reflects speciation (orthologs) from genes derived from duplication (paralogs). **HOGs** (Hierarchical Orthologous Groups; Altenhoff 2013) are the modern unit of analysis -- a maximal cluster of genes descended from a single ancestral gene at a defined taxonomic level. **OrthoFinder3** (Emms & Kelly 2025 bioRxiv) is the current Quest-for-Orthologs benchmark leader; its HOG output is 12% more accurate than v2 orthogroups (+20% with outgroup). For prokaryote-eukaryote scale, **TOGA** (Kirilenko 2023 Science) integrates whole-genome-alignment-based orthology with gene-loss classification.

The "ortholog conjecture" (orthologs more functionally similar than paralogs) is supported but weakly (Altenhoff 2012; Stamboulian 2020). Don't treat 1:1 ortholog labeling as automatic functional equivalence.

## Prerequisites

```bash
# Modern tools
conda install -c bioconda orthofinder=3.0
conda install -c bioconda sonicparanoid broccoli proteinortho

# OMA / FastOMA
wget https://omabrowser.org/standalone/OMA.tgz && tar xf OMA.tgz
pip install fastoma

# eggNOG-mapper
conda install -c bioconda eggnog-mapper

# JustOrthologs (close species, fast)
git clone https://github.com/ridgelab/justOrthologs

# TOGA (WGA-anchored orthology + gene loss)
conda env create -f https://raw.githubusercontent.com/hillerlab/TOGA/master/toga_env.yml

# Quality control
conda install -c bioconda busco compleasm
```

## Quick Start

Tell the AI agent what to do:
- "Run OrthoFinder3 on these proteomes with outgroup for HOG inference"
- "Cross-validate with SonicParanoid2 for consensus single-copy orthologs"
- "Identify ortholog vs paralog (in-paralog vs co-ortholog) for gene X across vertebrates"
- "Project annotations from human to chimp using TOGA with intactness classification"

## Example Prompts

### Standard HOG Inference

> "Run OrthoFinder3 on these 30 proteomes (vertebrates with outgroup invertebrate). Pre-filter to longest isoform per gene. Extract single-copy HOGs for phylogenomic concatenation. Cross-validate with SonicParanoid2 and report consensus single-copy set."

### Functional Annotation Transfer

> "Transfer functional annotations from human Ensembl release 113 to a newly sequenced primate genome. Use eggNOG-mapper for primary annotation + OrthoFinder3 HOGs for cross-validation. Apply ortholog-conjecture-aware confidence stratification: 1:1 ortholog + low dN/dS + EXP evidence = high confidence."

### Vertebrate-Scale Orthology with Gene Loss

> "Identify orthology and gene-loss across 100 mammal genomes using TOGA + CESAR. First build Cactus WGA. Then TOGA Nextflow pipeline. Report per-genome intactness distribution (I/PI/UL/L/M/PM). Identify mammalian lineages with notable gene loss."

## What the Agent Will Do

1. **Validate inputs**: BUSCO/Compleasm completeness per species; consistent annotation pipeline
2. **Pre-filter** to longest isoform per gene (AGAT or OrthoFinder primary_transcript.py)
3. **Include outgroup** at next-higher taxonomic level (essential for HOG rooting)
4. **Run primary tool**: OrthoFinder3 default for most cases; SonicParanoid2 for cross-validation
5. **Extract HOG at appropriate level**: N0 root for phylogenomic; per-clade levels for clade-specific
6. **Filter single-copy** for concatenation: present in >= 90% species, 1 copy each
7. **Classify orthology**: 1:1, 1:many, many:many; in-paralog vs out-paralog
8. **For functional transfer**: eggNOG-mapper + 1:1 + dN/dS < 0.2 + EXP evidence
9. **For WGA-anchored** (vertebrate-scale): TOGA + CESAR
10. **Report**: per-species classification, orthogroup statistics, consensus single-copy set
11. **Caveats**: WGD effects, ortholog conjecture, annotation heterogeneity, hidden paralogy

## Tips

- OrthoFinder3 (2025 bioRxiv) is the modern Quest-for-Orthologs leader
- HOG output (v3) is at `Phylogenetic_Hierarchical_Orthogroups/N0.tsv`; v2 layout differs
- Include >= 1 outgroup at next-higher taxonomic level for STRIDE rooting (+20% HOG accuracy)
- Pre-filter to longest isoform per gene; otherwise isoforms inflate co-ortholog counts
- For WGD-affected lineages, use synteny-aware methods (GENESPACE, ProteinOrtho-synteny)
- Cross-validate with 2+ methods (OrthoFinder3 + SonicParanoid2 or OrthoFinder3 + OMA)
- BUSCO/Compleasm completeness >= 90% required for inclusion
- TOGA (Kirilenko 2023) for vertebrate-scale orthology with explicit gene-loss classification
- Annotation pipeline normalization across species; mixing pipelines inflates accessory
- Ortholog conjecture is weakly supported; treat 1:1 ortholog functional transfer with care
- For HGT-affected lineages, OrthoFinder may return "vertical orthologs" that aren't real; use ALE
- eggNOG-mapper integrates orthology + functional annotation; current eggNOG database mandatory
- JustOrthologs is fastest for closely related species (same family); OrthoFinder3 for broader

## Related Skills

comparative-genomics/synteny-analysis - Synteny-anchored disambiguation in WGD lineages
comparative-genomics/whole-genome-duplication - Distinguishing ohnologs from orthologs
comparative-genomics/gene-family-evolution - CAFE5 birth-death on OG counts
comparative-genomics/gene-tree-species-tree-reconciliation - DTL-aware orthology for prokaryotes
comparative-genomics/positive-selection - Selection on single-copy alignments
comparative-genomics/comparative-annotation-projection - TOGA / CESAR uses orthology
phylogenetics/modern-tree-inference - Single-copy concatenation phylogenomics
phylogenetics/species-trees - Coalescent species tree from gene trees
alignment/multiple-alignment - MSA quality affects HOG output
genome-annotation/functional-annotation - eggNOG / orthology-based annotation
read-qc/quality-reports - BUSCO / Compleasm completeness

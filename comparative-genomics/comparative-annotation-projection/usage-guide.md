# Comparative Annotation Projection - Usage Guide

## Overview

Comparative annotation projection transfers gene annotations from a reference to query genome(s) using evolutionary conservation, replacing de novo prediction for many comparative-genomics workflows. The 2023-era standard is **TOGA + CESAR 2.0** (Kirilenko 2023 Science 380:eabn3107), which uses whole-genome alignment chains + ML classification + codon-aware exon projection to scale to hundreds of genomes (Zoonomia: 488 mammals; Bird10000 Genomes: 501 birds). **LiftOff** (Shumate & Salzberg 2021) is the standard for fast pairwise transfer without WGA.

The critical decision is **WGA-anchored** (TOGA: explicit gene-loss classification) vs **ortholog-anchored** (LiftOff: faster, more permissive).

## Prerequisites

```bash
# TOGA (WGA-based projection)
conda env create -f https://raw.githubusercontent.com/hillerlab/TOGA/master/toga_env.yml
# Requires Nextflow + nf-core for pipeline

# CESAR 2.0 (bundled with TOGA)
git clone https://github.com/hillerlab/CESAR2.0

# LiftOff (ortholog-based)
pip install liftoff

# UCSC liftOver (coordinate transfer)
conda install -c bioconda ucsc-liftover

# Comparative Annotation Toolkit (CAT)
git clone https://github.com/ComparativeGenomicsToolkit/Comparative-Annotation-Toolkit

# Cactus + HAL (for TOGA input)
docker pull quay.io/comparative-genomics-toolkit/cactus:latest

# De novo alternatives
conda install -c bioconda braker3 funannotate maker
```

## Quick Start

Tell the AI agent what to do:
- "Project mouse gene annotations onto rat using LiftOff"
- "Run TOGA + CESAR on a Cactus HAL to annotate 100 mammal genomes with gene-loss classification"
- "Build a comparative annotation pipeline with CAT (Snakemake) for multi-species annotation"
- "Identify pseudogenized genes in cetacean genomes using TOGA intactness classification"

## Example Prompts

### Vertebrate-Scale Annotation Projection

> "Annotate 100 mammal genomes using the human reference (GRCh38 + Ensembl 113). First build a Cactus WGA with the provided guide tree. Then run TOGA Nextflow pipeline to project the human annotation, producing per-gene intactness classification (I/PI/UL/L/M/PM). Identify genes lost in specific lineages (e.g. taste receptors in cetaceans)."

### Pairwise Annotation Transfer

> "Transfer Arabidopsis thaliana annotations to a newly sequenced Arabidopsis lyrata genome using LiftOff. Use `-copies` for tandem-rich regions (NLR clusters). Validate the projection by computing BUSCO on the resulting GFF."

### Multi-Reference Comparative Annotation

> "Build a comparative annotation pipeline for 20 plant genomes using CAT. Use Arabidopsis as primary reference + Solanum lycopersicum as secondary. Integrate de novo BRAKER3 predictions with LiftOff projections. Report per-species annotation quality metrics."

## What the Agent Will Do

1. **Validate inputs**: reference annotation BUSCO; assembly N50; species divergence < 300 Myr (TOGA limit)
2. **Build WGA** with Cactus (for TOGA) or skip (for LiftOff)
3. **Run projection** with appropriate tool:
   - TOGA: WGA-anchored with gene-loss classification
   - LiftOff: ortholog-anchored pairwise
   - CAT: integrated multi-reference pipeline
4. **Quality-check output**: BUSCO on projected annotation; per-gene intactness; pseudogenization
5. **Cross-validate** with RNA-Seq evidence for projected genes
6. **Document**: reference annotation version, divergence to query, projection success rate
7. **Caveats**: maximum divergence, tandem duplicate handling, reference bias, splice variant projection

## Tips

- TOGA is the Zoonomia / Bird10000 standard for clade-level annotation
- LiftOff is faster for pairwise transfer; doesn't require WGA
- Maximum divergence for TOGA: ~300 Myr (vertebrate); validate per clade
- Maximum divergence for LiftOff: ~80% nucleotide identity
- Tandem-rich regions (NLR clusters): LiftOff `-copies`; or pre-collapse for TOGA
- Reference annotation quality matters; verify BUSCO before propagating
- TOGA intactness codes (loss_summ_data.tsv): I (intact), PI (partial intact), UL (uncertain loss), L (lost), M (missing), PM (partial missing). Orthology classes (orthology_classification.tsv): one2one, one2many, many2one, many2many, PG (paralogous projection).
- For pseudogenization detection, combine TOGA + RNA-Seq + Ribo-Seq evidence
- For polyploid query, assign subgenomes first ([[whole-genome-duplication]])
- Multi-reference projection (CAT): consensus annotation; documented disagreements
- For fragmented assemblies, projection produces fragmented gene models; document
- Splice variants project independently; some isoforms may fail
- CESAR 2.0 is used internally by TOGA; standalone use rare
- For pairwise annotation transfer between closely related strains, LiftOff is gold standard
- GeMoMa integrates multiple reference species evidence; alternative for evidence-rich projection

## Related Skills

comparative-genomics/whole-genome-alignment - Cactus precedes TOGA
comparative-genomics/synteny-analysis - Synteny detection from WGA
comparative-genomics/ortholog-inference - TOGA orthology classification
comparative-genomics/pangenome-analysis - PGR-TK for repetitive / clinical genes
comparative-genomics/whole-genome-duplication - Subgenome assignment for polyploid query
genome-annotation/eukaryotic-gene-prediction - BRAKER3 / Funannotate de novo
genome-annotation/functional-annotation - Function assignment downstream
genome-annotation/annotation-transfer - Related skill on annotation transfer
read-qc/rnaseq-qc - RNA-Seq evidence to validate
read-alignment/star-alignment - RNA-Seq alignment for validation

# Eukaryotic Gene Prediction - Usage Guide

## Overview

Predict protein-coding gene structures in eukaryotic genomes by matching the pipeline to the available evidence: BRAKER3 (RNA-seq + protein), BRAKER1 (RNA-seq), GALBA or BRAKER2 (protein-only), Funannotate (fungi), GeMoMa (homology projection from a close relative), or Helixer/Tiberius (deep-learning ab initio when no evidence exists). The decisions that matter most are upstream of the tool: soft-mask repeats first, decontaminate and check assembly contiguity before training, and pick the right OrthoDB clade partition. Output is gene models in GFF3/GTF with protein and CDS sequences.

## Prerequisites

```bash
# BRAKER3 (wraps AUGUSTUS + GeneMark-ETP + TSEBRA) - prefer the official container for reproducibility
conda install -c bioconda braker3

# Protein-only and homology routes
conda install -c bioconda galba

# Fungal end-to-end pipeline
conda install -c bioconda funannotate

# Evidence prep + evaluation
conda install -c bioconda hisat2 star samtools busco

# Python parsing
pip install gffutils pandas
```

Requires a **soft-masked** genome assembly (run repeat-annotation first). Pin the OrthoDB partition version and the repeat library for reproducibility.

## Quick Start

Tell your AI agent what you want to do:
- "Predict genes in my eukaryotic genome using RNA-seq and protein evidence"
- "Run BRAKER3 on my soft-masked plant genome and check the gene-model sanity metrics"
- "I only have protein evidence from a close relative - which pipeline should I use?"
- "Annotate my fungal genome with Funannotate including UTRs and isoforms"

## Example Prompts

### Evidence-Matched Prediction

> "Run BRAKER3 on my soft-masked genome with these RNA-seq BAM files and the Vertebrata OrthoDB partition, and report mono-exonic fraction and protein-length distribution."

> "I have no RNA-seq and only distant proteins - run BRAKER2 with a broad OrthoDB clade, or Tiberius ab initio since this is a vertebrate."

### Protein-Only and Fungi

> "I have proteomes from two close relatives - run GALBA on my assembly."

> "Annotate my fungal genome with Funannotate, then run the update step to add UTRs and isoforms."

### Diagnosing a Poor Annotation

> "My annotation has 45,000 genes for a vertebrate - check whether this is haplotigs, unmasked repeats, or split models."

> "Compare assembly-BUSCO to proteome-BUSCO to tell whether my predictor missed present genes or the assembly is incomplete."

## What the Agent Will Do

1. Confirm the assembly is soft-masked (repeats in lowercase) and decontaminated, and check contiguity before training
2. Choose the pipeline from available evidence (BRAKER3 > BRAKER1/2 > GALBA/Funannotate/GeMoMa > DL ab initio)
3. Prepare evidence (splice-aware RNA-seq alignment; smallest OrthoDB partition containing the clade)
4. Run the predictor from its container with the correct species/clade parameters
5. Extract protein and CDS sequences
6. Run the triage panel: gene count vs relative, mono-exonic fraction, protein-length distribution, mRNA:gene ratio, BUSCO genome-vs-proteome
7. Flag failure signatures and recommend a PASA/Iso-Seq update step if isoforms/UTRs are needed

## Tips

- **Soft-masking is essential** - Hard-masking truncates real repeat-overlapping genes; soft-masking lets a real gene span a repeat. Filter the repeat library against a protein DB so NLR/zinc-finger families survive.
- **Training quality dominates** - A finder trained on a contaminated/fragmented assembly is confidently wrong genome-wide and BUSCO-green. Decontaminate and check assembly BUSCO/N50 first.
- **OrthoDB clade** - Pick the smallest partition that still contains the clade (Vertebrata for a fish, not all Metazoa). GALBA wants a few close-relative proteomes instead.
- **One isoform, no UTRs by default** - If `mRNA:gene == 1.00`, the annotation is isoform/UTR-naive; alternative-splicing and 3'-tag scRNA-seq analyses against it are untrustworthy. Add a PASA/Iso-Seq update step.
- **Fungi -> Funannotate** - It is the de facto fungal standard (tiny-intron-aware, bundled EVM + PASA update); BRAKER's `--fungus` flag exists but Funannotate is the smoother fungal path.
- **High BUSCO-Duplicated** - >5-8% with no known WGD means purge_dups the assembly first, not "deduplicate genes" afterward. In a known polyploid, high D is real - confirm via synteny/Ks and keep.
- **Reproducibility** - Run from official containers; the GFF3 alone is not reproducible without the masked genome, evidence, trained parameters, and versions.

## Related Skills

- genome-annotation/repeat-annotation - Soft-mask repeats before gene prediction
- genome-annotation/functional-annotation - Add functional annotation to predicted genes
- genome-annotation/annotation-qc - BUSCO genome-vs-proteome, OMArk, gene-set sanity
- read-alignment/star-alignment - Splice-aware RNA-seq aligner for evidence
- genome-assembly/assembly-qc - Verify assembly quality and purge haplotigs first

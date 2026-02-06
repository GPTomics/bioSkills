# Repeat Annotation - Usage Guide

## Overview

Identify, classify, and mask repetitive elements and transposable elements in genome assemblies using RepeatModeler (de novo repeat library construction) and RepeatMasker (genome-wide repeat annotation and masking). Softmasked output from RepeatMasker is a prerequisite for eukaryotic gene prediction with BRAKER3.

## Prerequisites

```bash
# RepeatModeler and RepeatMasker
conda install -c bioconda repeatmodeler repeatmasker

# TEtranscripts for TE expression quantification
pip install TEtranscripts

# Python utilities
pip install pandas matplotlib biopython
```

## Quick Start

Tell your AI agent what you want to do:
- "Mask repeats in my genome assembly before gene prediction"
- "Build a de novo repeat library and annotate transposable elements"
- "Quantify transposable element expression from RNA-seq"

## Example Prompts

### Repeat Masking

> "Run RepeatModeler to build a de novo repeat library, then RepeatMasker to softmask my genome"

> "Mask repeats in my assembly using the human Dfam library"

### TE Analysis

> "Classify all transposable elements in my genome and report percentages"

> "Generate a repeat divergence landscape plot"

### TE Expression

> "Quantify TE expression in treated vs control samples using TEtranscripts"

> "Run differential TE expression analysis alongside gene expression"

## What the Agent Will Do

1. Build RepeatModeler database from the assembly
2. Run RepeatModeler to construct de novo repeat consensus library
3. Run RepeatMasker with the de novo library to annotate and softmask repeats
4. Generate repeat statistics (content by class, divergence distribution)
5. Produce softmasked FASTA ready for gene prediction
6. Optionally quantify TE expression with TEtranscripts

## Tips

- **Always run RepeatModeler first** - De novo libraries capture species-specific repeats missed by curated databases
- **Combine libraries** - Use both de novo and Dfam/RepBase libraries for best results
- **Softmasking** - Always use `-xsmall` flag; hardmasking (N's) loses sequence information
- **Gene prediction** - Softmasked genome is required input for BRAKER3 and most gene predictors
- **STAR alignment for TEtranscripts** - Must retain multi-mapping reads (--outFilterMultimapNmax 100)
- **Large genomes** - RepeatModeler takes days for mammalian-size genomes; plan accordingly
- **EDTA** - Consider EDTA as alternative to RepeatModeler for plant genomes with many LTR elements

## Related Skills

- genome-annotation/eukaryotic-gene-prediction - Requires softmasked genome
- genome-assembly/assembly-qc - Assembly quality including repeat content
- differential-expression/deseq2-basics - Differential TE expression

# Ortholog Inference - Usage Guide

## Overview

Infer orthologous gene groups across species using OrthoFinder and ProteinOrtho for comparative genomics and functional annotation transfer.

## Prerequisites

```bash
# OrthoFinder (recommended)
conda install -c bioconda orthofinder

# ProteinOrtho (faster alternative)
conda install -c bioconda proteinortho

# DIAMOND (used by both)
conda install -c bioconda diamond
```

## Quick Start

Tell your AI agent what you want to do:
- "Find orthologs of BRCA1 across vertebrates"
- "Build orthogroups from these proteomes"
- "Extract single-copy orthologs for phylogenomics"

## Example Prompts

### Orthogroup Analysis

> "Run OrthoFinder on my proteome files"

> "How many single-copy orthologs are there across all species?"

### Specific Gene Queries

> "Find orthologs of this gene in mouse and zebrafish"

> "Identify paralogs within each species"

### Phylogenomics

> "Extract single-copy orthologs for building a species tree"

> "Get universal orthologs present in all genomes"

## What the Agent Will Do

1. Prepare proteome FASTA files (one per species)
2. Run OrthoFinder or ProteinOrtho for ortholog detection
3. Parse orthogroup results into structured format
4. Classify orthogroups (single-copy, universal, partial)
5. Extract single-copy orthologs if needed
6. Identify paralogs and co-orthologs
7. Enable annotation transfer via orthology

## Tips

- **Input format** - One proteome FASTA per species; filename becomes species name
- **Single-copy** - Ideal for phylogenomics; one gene per species, no paralogs
- **OrthoFinder** - More accurate gene trees; better for evolution studies
- **ProteinOrtho** - Faster; good for many genomes or quick surveys
- **Co-orthologs** - Multiple genes orthologous to one gene indicate duplication
- **Annotation transfer** - Most reliable for single-copy orthologs

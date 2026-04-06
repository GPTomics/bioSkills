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

1. Verify proteome completeness with BUSCO/Compleasm
2. Prepare proteome FASTA files (one per species, isoforms removed)
3. Select method based on dataset size and accuracy needs (OrthoFinder vs ProteinOrtho)
4. Run ortholog detection with appropriate parameters
5. Parse orthogroup results and classify by copy number pattern
6. Extract single-copy orthologs for phylogenomics if needed
7. Identify in-paralogs, out-paralogs, and co-orthologs
8. Enable annotation transfer with appropriate confidence levels

## Tips

- **Input quality** - Remove isoforms (keep longest per gene) to avoid inflating copy numbers; verify completeness with BUSCO first
- **Annotation consistency** - Heterogeneous annotations across species create false lineage-specific expansions; use consistent pipelines when possible
- **OrthoFinder** - Tree-based, most accurate; use `-M msa` for <20 species; default for evolutionary analysis
- **ProteinOrtho** - Graph-based, faster; good for 50+ genomes or quick surveys
- **OMA/FastOMA** - Highest precision but lowest recall; use when false positives are costly
- **Single-copy orthologs** - Ideal for phylogenomics; one gene per species, no paralogy complications
- **In-paralogs vs out-paralogs** - Distinguishing them requires speciation context; OrthoFinder resolves via gene tree reconciliation
- **Annotation transfer** - Highest confidence for 1:1 orthologs; decreases with co-orthologs and distant homologs

## Related Skills

- comparative-genomics/synteny-analysis - Synteny-based ortholog verification
- comparative-genomics/positive-selection - Selection analysis on orthologs
- phylogenetics/modern-tree-inference - Build trees from single-copy orthologs
- genome-annotation/annotation-transfer - Transfer annotations via orthology

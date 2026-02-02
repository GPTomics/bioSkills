# HGT Detection - Usage Guide

## Overview

Detect horizontal gene transfer events using HGTector, compositional analysis, and phylogenetic incongruence methods for prokaryotic genome evolution studies.

## Prerequisites

```bash
# HGTector
pip install hgtector

# Reference databases
hgtector database  # Download and setup

# IQ-TREE for topology tests
conda install -c bioconda iqtree

# Python dependencies
pip install biopython pandas numpy
```

## Quick Start

Tell your AI agent what you want to do:
- "Detect horizontally transferred genes in this bacterial genome"
- "Find genomic islands with anomalous GC content"
- "Identify genes with unexpected phylogenetic placement"

## Example Prompts

### Compositional Analysis

> "Find genes with anomalous GC content in this genome"

> "Calculate codon usage bias for potential foreign genes"

### HGTector Analysis

> "Run HGTector to find putative HGT events"

> "Which genes have unusual taxonomic distribution?"

### Phylogenetic Methods

> "Test for phylogenetic incongruence in this gene family"

> "Does this gene tree conflict with the species tree?"

### Genomic Islands

> "Identify genomic islands in this bacterial genome"

> "Find clusters of foreign genes with mobile elements"

## What the Agent Will Do

1. Calculate genome-wide GC content baseline
2. Identify genes with anomalous composition (GC, codon usage)
3. Run HGTector for phyletic distribution analysis
4. Flag genes with unexpected taxonomic matches
5. Cluster anomalous genes into genomic islands
6. Annotate islands for mobile elements
7. Report HGT candidates with confidence levels

## Tips

- **GC threshold** - Z-score > 2 suggests foreign origin; adjust for genome
- **Codon usage** - CAI < 0.5 indicates potential HGT
- **Island size** - Minimum 3 genes for confident island calls
- **Mobile elements** - Integrases/transposases support HGT origin
- **tRNA sites** - Common integration hotspots for genomic islands
- **Amelioration** - Ancient HGT may have host-like composition

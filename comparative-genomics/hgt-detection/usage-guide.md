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

1. Calculate genome-wide GC content and codon usage baselines
2. Identify genes with anomalous composition (GC z-score, CAI)
3. Run HGTector for phyletic distribution analysis if database available
4. Flag genes with unexpected taxonomic distribution
5. Cluster anomalous genes into genomic islands
6. Annotate islands for mobile element signatures (integrases, transposases)
7. Cross-validate with phylogenetic incongruence where possible
8. Report HGT candidates with multi-method confidence levels (strong/moderate/suggestive)

## Tips

- **Use multiple methods** - No single method catches all HGT; combine compositional + phylogenetic + HGTector for robust calls
- **GC threshold** - |Z| > 2 moderate, |Z| > 3 strong; adjust for genomes with high GC variance (e.g., Streptomyces)
- **Amelioration window** - Compositional methods only detect recent HGT (<10 Myr); ancient transfers are fully ameliorated
- **ILS vs HGT** - Gene tree discordance can be incomplete lineage sorting, not HGT; expect ILS at rapid radiations
- **Contamination** - Always rule out assembly contamination before claiming HGT, especially in draft assemblies
- **Island signatures** - Minimum 3 genes, flanking mobile elements, tRNA integration sites all increase confidence
- **Eukaryotic HGT** - Rarer but documented; exclude organellar gene transfer (EGT) as alternative
- **Donor identification** - Phylogenetic placement within a distant clade (not sister) is the strongest HGT signal

## Related Skills

- comparative-genomics/ortholog-inference - Identify orthologs for phylogenetic tests
- phylogenetics/modern-tree-inference - Build gene trees for incongruence analysis
- metagenomics/amr-detection - AMR genes often on mobile elements

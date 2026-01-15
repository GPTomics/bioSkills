# Sequence Properties - Usage Guide

## Overview

This skill enables AI agents to help you calculate physical and chemical properties of biological sequences using Biopython. It covers GC content, molecular weight, melting temperature, and protein properties like isoelectric point and hydropathicity.

## Prerequisites

```bash
pip install biopython
```

## Quick Start

Tell your AI agent what you want to do:

- "Calculate the GC content of this sequence"
- "What is the molecular weight of this DNA?"
- "Calculate the melting temperature for this primer"
- "Analyze this protein sequence for stability and pI"

## Example Prompts

### GC Content
> "What is the GC content of each sequence in my FASTA file?"

### Molecular Weight
> "Calculate the molecular weight of this double-stranded DNA"

### Melting Temperature
> "What's the Tm of this 20-mer primer?"

### Protein Analysis
> "Give me a full analysis of this protein: MW, pI, instability index, GRAVY"

### Sliding Window
> "Plot the GC content along this sequence using 100 bp windows"

### Codon Usage
> "What is the codon usage bias in this coding sequence?"

## What the Agent Will Do

1. Import appropriate modules (SeqUtils, ProtParam)
2. Parse the sequence
3. Calculate requested properties
4. Return formatted results
5. For multiple sequences, summarize statistics

## DNA/RNA Properties Available

- **GC content**: Fraction of G+C bases
- **Molecular weight**: Mass in Daltons
- **Melting temperature**: Tm for hybridization
- **Base composition**: Percentage of each nucleotide
- **Dinucleotide frequencies**: CpG, etc.

## Protein Properties Available

- **Molecular weight**: Mass in Daltons
- **Isoelectric point (pI)**: pH of zero net charge
- **Instability index**: Stability prediction (<40 = stable)
- **GRAVY**: Hydropathicity score
- **Aromaticity**: Fraction of aromatic amino acids
- **Secondary structure**: Helix/sheet/turn fractions
- **Extinction coefficient**: For concentration measurements

## Tips

- GC content returns a fraction (0-1), multiply by 100 for percentage
- For accurate Tm, use the nearest-neighbor method (Tm_NN)
- Protein analysis requires standard amino acid letters only
- For genome-wide analysis, use sliding windows

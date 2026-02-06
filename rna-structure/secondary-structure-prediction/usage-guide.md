# Secondary Structure Prediction - Usage Guide

## Overview
Predict RNA secondary structures from sequence using thermodynamic models. ViennaRNA computes minimum free energy (MFE) structures, partition function ensembles, base-pair probabilities, consensus structures from alignments, and RNA-RNA interaction structures.

## Prerequisites
```bash
# ViennaRNA (includes RNAfold, RNAalifold, RNAcofold, Python API)
conda install -c bioconda viennarna

# Optional: LinearFold for long sequences
conda install -c bioconda linearfold

# Python dependencies
pip install matplotlib numpy
```

## Quick Start
Tell your AI agent what you want to do:
- "Predict the secondary structure of this RNA sequence"
- "Fold this tRNA sequence and show the base-pair probabilities"
- "Build a consensus structure from my alignment of homologous RNAs"
- "Compare MFE structures between wild-type and mutant sequences"
- "Use SHAPE reactivities to constrain RNA folding"

## Example Prompts

### Single Sequence Folding
> "Fold this RNA sequence and report the MFE, centroid, and MEA structures."

> "Predict the structure of my 5' UTR and tell me if it has a stable hairpin."

### Ensemble Analysis
> "Compute the partition function for my RNA and show the base-pair probability dot plot."

> "How well-defined is the structure of this RNA? Show the ensemble diversity."

### Consensus Structures
> "I have a Stockholm alignment of my RNA family. Predict the consensus structure with covariation."

> "Compare thermodynamic folding with covariation-supported base pairs in my alignment."

### RNA-RNA Interaction
> "Predict the hybridization structure between my sRNA and its target mRNA."

> "What is the binding energy between these two RNA sequences?"

### Constrained Folding
> "Fold my RNA but force positions 15-20 to be unpaired."

> "Use my SHAPE reactivity data to guide secondary structure prediction."

## What the Agent Will Do
1. Select appropriate ViennaRNA tool based on input type (single sequence, alignment, two sequences)
2. Run folding with recommended parameters (partition function, no lonely pairs if appropriate)
3. Parse output for MFE, structure, and base-pair probabilities
4. Visualize results (dot plot, arc diagram, or forna-compatible JSON)
5. Assess structure confidence using ensemble diversity and base-pair probabilities

## Tips
- **Partition function** - Always use `-p` to get base-pair probabilities alongside MFE; ensemble information is more informative than MFE alone
- **Long sequences** - For RNAs longer than 5,000 nt, use LinearFold or `--maxBPspan` to limit folding window
- **Pseudoknots** - ViennaRNA does not predict pseudoknots; use pKiss or IPknot if pseudoknots are expected
- **Temperature** - Default 37C is correct for most applications; adjust with `-T` for in vitro experiments at different temperatures
- **Lonely pairs** - Use `--noLP` to suppress isolated base pairs (improves prediction for structured RNAs)
- **Covariation** - RNAalifold consensus structures are more reliable than single-sequence predictions for conserved RNAs

## Related Skills
- ncrna-search - Classify structured RNAs by family using Infernal/Rfam
- structure-probing - Experimental SHAPE/DMS data to constrain predictions
- genome-annotation/ncrna-annotation - Genome-wide ncRNA annotation
- sequence-manipulation/sequence-properties - Sequence composition and GC content

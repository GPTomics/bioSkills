# Secondary Structure Prediction - Usage Guide

## Overview
Predict RNA secondary structure with ViennaRNA, reporting the Boltzmann ensemble (partition function, base-pair probabilities, centroid, MEA, per-base confidence, stochastic samples) rather than a single minimum free energy (MFE) fold. Covers consensus folding from alignments, SHAPE-constrained folding, RNA-RNA interaction, local and linear-time methods for long RNA, and the cases where pseudoknot-aware or deep-learning tools are needed.

## Prerequisites
```bash
# ViennaRNA (RNAfold, RNAalifold, RNAcofold, RNAduplex, RNAup, RNAplfold, RNALfold, Python API)
conda install -c bioconda viennarna

# Optional: linear-time folding for long sequences
conda install -c bioconda linearfold linearpartition

# Optional: pseudoknot-aware folding
conda install -c bioconda ipknot

# Python dependencies
pip install matplotlib numpy
```

## Quick Start
Tell your AI agent what you want to do:
- "Fold this RNA and report the ensemble, not just the MFE"
- "How well-defined is this structure? Show me the per-base confidence"
- "Build a consensus structure from my alignment of homologs"
- "This RNA is 8 kb, fold it sensibly"
- "My RNA probably has a pseudoknot, predict it"
- "Use my SHAPE reactivities to constrain folding"

## Example Prompts

### Single Sequence Folding
> "Fold this RNA and report the MFE, centroid, and MEA structures, and tell me which to trust."

> "Predict the structure of my 5' UTR and tell me whether the hairpin is well-defined or ambiguous."

### Ensemble Confidence
> "Compute the base-pair probability dot plot and the positional entropy for this RNA."

> "Is this RNA likely to adopt a single structure or switch between conformations?"

### Consensus Structures
> "I have a Stockholm alignment of my RNA family. Predict the consensus structure with covariation, and tell me whether the covariation is statistically significant."

### Long RNA
> "Fold this 10 kb mRNA and give me local accessibility, not one global structure."

### RNA-RNA Interaction
> "Predict the binding between my sRNA and its target mRNA, accounting for whether the target site is accessible."

### Constrained / Pseudoknot / SHAPE
> "Fold my RNA but force positions 15-20 unpaired."

> "My RNA is a frameshift element with a pseudoknot; standard folding misses it."

> "Use my SHAPE reactivity data to guide secondary structure prediction."

## What the Agent Will Do
1. Choose the ViennaRNA program by the question (single sequence, alignment, two strands, long RNA, pseudoknot)
2. Run partition-function folding (`-p`) so ensemble quantities are available, not just the MFE
3. Report the appropriate answer (MFE/centroid/MEA/sampling) with confidence (base-pair probabilities, positional entropy, ensemble diversity)
4. Steer long RNAs to local/linear methods and suspected pseudoknots to pseudoknot-aware tools
5. Flag when a predicted structure needs covariation (R-scape) or probing validation

## Tips
- **Ensemble first** - Always run `-p`; the base-pair probability matrix and positional entropy say which parts of the fold are trustworthy. The MFE alone hides the uncertainty and is often not the functional structure.
- **Pick the answer deliberately** - Centroid is conservative (fewer false pairs); MEA tunes precision/recall via gamma; stochastic sampling represents alternative conformations. Report MFE alone only as a quick first look.
- **Long RNA** - A single global MFE for a multi-kilobase RNA is nearly meaningless; use RNAplfold (local/accessibility), LinearFold, or LinearPartition. Folding degrades sharply past ~700 nt.
- **Pseudoknots** - ViennaRNA cannot predict pseudoknots (it returns the best nested fold). For frameshift elements, riboswitch aptamers, telomerase, and viral UTRs, use IPknot, ProbKnot, or Knotty.
- **Comparative is a hypothesis** - RNAalifold consensus structures beat single-sequence folding only when the alignment carries real covariation; a predicted conserved structure is not established until R-scape shows significant covariation (HOTAIR/Xist/SRA failed this test).
- **Deep learning is not a default** - End-to-end DL predictors do not reliably generalize to novel families; for unseen RNAs prefer thermodynamics plus covariation/probing, or the hybrid MXfold2.
- **PostScript clutter** - RNAfold/RNAalifold write `*_dp.ps`/`*_ss.ps` to the working directory; pass `--noPS`.

## Related Skills
- structure-probing - Obtain SHAPE/DMS reactivities to constrain folding
- ncrna-search - Classify structured RNAs by family using Infernal/Rfam
- covariation-analysis - Statistically validate a predicted conserved structure with R-scape
- genome-annotation/ncrna-annotation - Genome-wide ncRNA annotation
- small-rna-seq/target-prediction - miRNA-target prediction using accessibility
- sequence-manipulation/sequence-properties - Sequence composition and GC content

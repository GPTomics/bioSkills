# Alignment Trimming - Usage Guide

## Overview

This skill covers post-alignment column filtering and per-residue cleaning of multiple sequence alignments. It guides tool selection (ClipKIT, trimAl, BMGE, Divvier, HMMcleaner, Gblocks) by downstream goal (phylogenetic-tree input, HMM profile building, selection analysis, structure modelling) rather than treating "trimming" as a single decision.

## Prerequisites

```bash
# CLI tools (all primary trimmers)
conda install -c bioconda clipkit trimal
# OR pip:
pip install clipkit
# trimAl: https://github.com/inab/trimal

# BMGE (Java)
# https://gitlab.pasteur.fr/GIPhy/BMGE

# Divvier (compile from source)
# https://github.com/simonwhelan/Divvier

# HMMcleaner (Perl + HMMER required)
# cpanm Bio::MUST::Apps::HmmCleaner
```

## Quick Start

Tell your AI agent what you want to do:
- "Trim this alignment for phylogenetic-tree input using ClipKIT"
- "Prepare an HMM-ready alignment with aggressive gap removal"
- "Remove cross-contaminating residues with HMMcleaner before phylogenomic analysis"
- "Generate a structural-MSA-friendly trimming that preserves codon boundaries"

## Example Prompts

### Phylogenetic Input Trimming
> "I have a concatenated supermatrix of 500 genes; what trimmer should I use for IQ-TREE input?"

> "Trim this single-gene alignment with ClipKIT smart-gap and report what fraction of columns were retained"

> "Compare ClipKIT kpic-smart-gap vs trimAl -strictplus on the same alignment and tell me which produces a more stable ML tree"

### HMM Profile Preparation
> "Trim this alignment for HMMER hmmbuild input using trimAl gappyout mode"

> "Strip insert columns and lowercase characters before building a profile HMM"

### Cross-Contamination Cleaning
> "I suspect some sequences have contaminating residues; clean them with HMMcleaner first then trim"

> "Mask low-confidence residues with HMMcleaner before building a tree"

### Selection-Analysis-Safe Filtering
> "Mark unreliable columns for codeml input but DO NOT remove them; selection analysis is sensitive to column removal"

> "Use TCS column scores to mask unreliable codon positions before BUSTED/MEME/aBSREL"

### Mode Comparison
> "Run all three ClipKIT modes (smart-gap, kpic-smart-gap, gappy) and compare retention fractions"

> "Test how aggressive each trimAl mode is on my alignment and recommend one"

## What the Agent Will Do

1. Identify the downstream goal (tree, HMM, selection analysis, motif scan, structure modelling)
2. Select the appropriate trimmer and mode (ClipKIT modes vs trimAl modes vs BMGE entropy threshold)
3. Run with reproducibility-preserving flags (column mapping, log file)
4. Compute trimming fraction and warn if > 20-30% of columns are removed
5. For selection analysis, prefer per-column reliability scores (TCS, GUIDANCE2) over column removal
6. Optionally compare alternative trimming modes on the same input

## Key Concepts

| Concept | Meaning |
|---------|---------|
| Column retention | Fraction of original columns kept; > 70% is usually safe, < 50% concerning |
| Aggressiveness | How readily a mode removes columns; ClipKIT smart-gap is moderate, trimAl strictplus is high |
| Column mapping | Index list of which original columns survived (`-colnumbering`, `--log`) |
| Per-residue masking | Replace bad residues with `X` instead of removing the column (HMMcleaner) |
| Column splitting | Divvier alternative: split ambiguous columns rather than remove them |
| Codon-aware trimming | Preserve frame boundaries for selection analysis (MACSE trimAlignment) |

## Tips

- ClipKIT `kpic-smart-gap` is the recommended default for phylogenetic-tree input on concatenated supermatrices; `smart-gap` for single genes
- trimAl `-gappyout` is the right tool for HMM profile preparation; aggressive gap removal benefits HMMER `hmmbuild`
- For selection analysis (PAML codeml, HyPhy), DO NOT aggressively trim; column removal causes false-positive selection signals (Fletcher & Yang 2010 MBE)
- BMGE is the prokaryotic phylogenomics standard (default in GToTree); recommended `-h 0.4 -g 0.2`
- Divvier preserves more phylogenetic signal than column removal because ambiguous columns are SPLIT rather than dropped
- HMMcleaner masks contaminating residues with `X`; pair with a column trimmer for double cleaning
- Always retain column-mapping (`-colnumbering`, `--log`) for reproducibility
- The 20%/40% rule reconciles the Tan 2015 vs Steenwyk 2020 controversy: <20% column removal is neutral; >40% removes phylogenetic signal alongside noise. If your trimmer drops more than 40% of columns, the mode is too aggressive for the dataset
- For phylogenomic incongruence due to suspected alignment artefact (rather than ILS or introgression), apply PhyIN (Maddison 2024 PeerJ) as a second-pass trimmer after ClipKIT/trimAl; it flags neighbouring columns whose split patterns are pairwise tree-incompatible
- Gblocks default parameters are too aggressive for modern use; relax them or switch to ClipKIT/trimAl
- Codon-aware trimming (MACSE `trimAlignment` or per-codon-block ClipKIT) is required for any post-MSA cleaning of dN/dS input

## Related Skills

- alignment/multiple-alignment - Generate the input MSA; confidence assessment (GUIDANCE2, TCS) for selection input
- alignment/msa-parsing - Parse and post-process trimmed alignments
- alignment/msa-statistics - Conservation and gap statistics to inform trimmer choice
- alignment/structural-alignment - Apply same trimming logic to structural MSAs
- phylogenetics/modern-tree-inference - Build trees from trimmed alignments

# Covalent Inhibitor Design Usage Guide

## Overview

Design covalent inhibitors targeting Cys, Lys, Ser, Thr, Tyr, or Asp residues. Balance warhead reactivity, GSH stability, geometric accessibility, and irreversible vs reversible covalent. Covers DOCKovalent, HCovDock, GOLD covalent, and reactivity-aware SAR.

## Prerequisites

```bash
pip install rdkit
# DOCKovalent: web service (covalent.docking.org)
# GOLD: commercial license
# HCovDock: standalone installation
```

## Quick Start

Tell the AI agent what to do:
- "Suggest warheads for cysteine-selective covalent inhibitor"
- "Score acrylamide-containing analogs for reactivity"
- "Plan reactive group SAR for KRAS G12C inhibitor"
- "Filter library for GSH-stable warheads (acrylamide alpha-substituted)"

## Example Prompts

### Cysteine targeting
> "Identify acrylamide and chloroacetamide candidates in library.smi. Filter by GSH stability surrogate (alpha-substituents on warhead). Output sorted by reactivity tier."

### Reactivity SAR
> "For 30 acrylamide compounds in series.csv, compute LUMO surrogate (alpha-C substituent count). Correlate with measured kinact/Ki."

### Covalent docking
> "Dock 50 acrylamide candidates against EGFR C797. Use HCovDock; require covalent geometry to Cys797 within 4 A. Rank by pose score + reactivity."

### Bivalent / PROTAC
> "Combine cysteine-warhead inhibitor with E3 ligand via 12-atom linker. Predict ternary complex stability."

## What the Agent Will Do

1. Identify warheads in input (acrylamide, chloroacetamide, vinyl sulfone, etc.).
2. Score reactivity (LUMO surrogate, alpha-substituents).
3. Run covalent docking if requested (DOCKovalent / HCovDock / GOLD).
4. Validate geometry: ligand warhead within reach of target Cys/Lys/Tyr.
5. Filter by GSH stability surrogate (off-target risk).
6. Output ranked candidates with covalent pose + reactivity tier.

## Tips

- Cysteine targets cover ~98% of covalent drugs.
- Acrylamide is the drug-development standard; chloroacetamide is too reactive for clinical.
- GSH t1/2 > 4 hours is the modern target for drug-like covalent.
- Always validate geometric fit: warhead must be within 4 A of nucleophilic residue.
- Use kinact/Ki not IC50 for ranking; covalent kinetics matter.

## Related Skills

- chemoinformatics/substructure-search - Warhead SMARTS detection
- chemoinformatics/virtual-screening - Non-covalent docking first
- chemoinformatics/pose-validation - Validate covalent docking
- chemoinformatics/molecular-descriptors - Reactivity surrogates
- chemoinformatics/admet-prediction - ADMET of covalent leads
- chemoinformatics/protac-degraders - PROTAC with covalent warhead

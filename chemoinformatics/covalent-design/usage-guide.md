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
- "Classify candidate warheads, then identify which compounds require experimental GSH-reactivity testing"

## Example Prompts

### Cysteine targeting
> "Identify acrylamide and chloroacetamide candidates in library.smi. Report alpha substitution as a structural feature and flag compounds for matched GSH-reactivity measurements."

### Reactivity SAR
> "For 30 acrylamide compounds in series.csv, compute alpha-C substituent count and compare it with measured GSH rates and kinact/Ki without treating it as a LUMO calculation."

### Covalent docking
> "Dock acrylamide candidates against EGFR C797 with a reaction-appropriate protocol. Inspect Cys797 Sγ-to-electrophile distance and approach geometry, then integrate measured reactivity."

### Bivalent / PROTAC
> "Combine cysteine-warhead inhibitor with E3 ligand via 12-atom linker. Predict ternary complex stability."

## What the Agent Will Do

1. Identify warheads in input (acrylamide, chloroacetamide, vinyl sulfone, etc.).
2. Compute structural features and integrate measured or validated reactivity data.
3. Run covalent docking if requested (DOCKovalent / HCovDock / GOLD).
4. Validate reaction-atom geometry: for cysteine, use Sγ and the ligand electrophilic atom, not Cβ.
5. Review experimental GSH reactivity and off-target evidence.
6. Output ranked candidates with covalent pose + reactivity tier.

## Tips

- Cysteine is a common covalent target, but prevalence depends on the drug set and counting method.
- Acrylamide and haloacetamide behavior depends on the complete compound and target context; do not assign suitability from warhead class alone.
- Report GSH assay conditions and measured kinetics rather than applying a universal half-life cutoff.
- Calibrate geometric criteria for the reaction and docking protocol; no universal distance cutoff replaces reaction-specific angles and atom identities.
- Use `kinact/Ki` for irreversible two-step systems that fit that model; use mechanism-appropriate reversible-covalent measurements otherwise.

## Related Skills

- chemoinformatics/substructure-search - Warhead SMARTS detection
- chemoinformatics/virtual-screening - Non-covalent docking first
- chemoinformatics/pose-validation - Validate covalent docking
- chemoinformatics/molecular-descriptors - Reactivity surrogates
- chemoinformatics/admet-prediction - ADMET of covalent leads
- chemoinformatics/protac-degraders - PROTAC with covalent warhead

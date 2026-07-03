# Structure Validation - Usage Guide

## Overview

This skill decides whether a macromolecular model, or a specific region of it, is reliable enough to build on before any downstream measurement, docking, or mechanistic claim. It treats a coordinate file as an interpreted model whose reliability varies per-atom, so it reads the local signals (per-residue B-factor, real-space fit) rather than trusting the single headline resolution. It covers experimental validation (resolution, R-work vs R-free and the overfitting gap, B-factor sanity, MolProbity clashscore/Ramachandran/rotamer outliers, cis non-proline peptides, the wwPDB report and its percentile sliders) and predicted-model validation (pLDDT bands, PAE inter-domain confidence, and trimming with phenix.process_predicted_model before molecular replacement or docking). It also handles the cryo-EM global-vs-local resolution distinction (FSC 0.143 half-map vs 0.5 map-model) and NMR ensemble spread.

## Prerequisites

```bash
pip install biopython numpy
```

MolProbity, phenix (`phenix.molprobity`, `phenix.process_predicted_model`), and DSSP (`mkdssp`) are separate CLI installs, not pip packages; install them independently and confirm they are on PATH. The MolProbity web service (http://molprobity.biochem.duke.edu) is an alternative to a local phenix install.

## Quick Start

Tell the AI agent what you want to decide:
- "Is this 2.8 Angstrom structure good enough to trust the active-site geometry?"
- "Read the resolution and R-free from this mmCIF and tell me if it is overfit"
- "Flag the Ramachandran and rotamer outliers and any cis non-proline peptides"
- "Which residues in this chain have suspiciously high B-factors?"
- "Validate this AlphaFold model before I dock a ligand into it"
- "Should I trust the domain arrangement in this predicted structure?"
- "Is the resolution uniform across this cryo-EM map or is the periphery a guess?"
- "Report the per-residue spread of this NMR ensemble instead of using model 1"

## Example Prompts

### Judging a deposited experimental structure
> "I want to measure a catalytic distance in this 3.1 Angstrom crystal structure. Read the resolution, R-work, R-free and the overfitting gap, then check whether the active-site residues have high local B-factors before I trust the measurement."

> "Run the geometry validation on this model and tell me the clashscore, Ramachandran favored percentage and outliers, and poor-rotamer percentage. Is this worse than typical for its resolution?"

> "Pull the wwPDB validation report for this PDB ID and tell me whether my region of interest is flagged as an RSRZ outlier."

### Validating a predicted model
> "Before I use this AlphaFold model for molecular replacement, trim the low-pLDDT residues, convert pLDDT to a pseudo-B, and split it into PAE-defined domains."

> "This predicted structure looks confident everywhere. Read the PAE matrix and tell me whether the two domains have a trustworthy relative orientation or are independently placed."

> "A long stretch of this model is low pLDDT. Is that a modeling failure or an intrinsically disordered region I should keep?"

### Cryo-EM and NMR
> "This is a 2.6 Angstrom cryo-EM structure of a large complex. Is that resolution uniform, or is the peripheral arm at much lower local resolution and effectively a docked hypothesis?"

> "Compute the per-residue Ca RMSD across all models of this NMR ensemble and show me which regions are well-restrained versus flexible. Do not average the coordinates."

## What the Agent Will Do

1. Identify the experimental method (`_exptl.method`) to route to the correct validation logic (X-ray, cryo-EM, NMR, or predicted).
2. Read resolution, R-work, and R-free from the mmCIF header via MMCIF2Dict and compute the R-free-minus-R-work overfitting gap against a resolution-scaled expectation.
3. Sanity-check per-residue B-factors within the structure and flag the least-certain residues (never comparing B across structures without normalization).
4. Screen backbone geometry for Ramachandran outliers and cis non-proline peptides in Python, then run phenix.molprobity for archive-calibrated clashscore, rotamer, and Ramachandran percentages when a real answer is needed.
5. For predicted models, read pLDDT into bands, read the PAE matrix for domain segmentation, and process with phenix.process_predicted_model before docking or molecular replacement.
6. For cryo-EM, direct the reader to the EMDB local-resolution map and distinguish the 0.143 (half-map) and 0.5 (map-model) FSC thresholds.
7. For NMR, compute and report the per-residue ensemble spread instead of averaging or picking model 1.

## Tips

- The single global resolution is a data ceiling, not a per-region quality certificate. A 1.5 Angstrom structure can still have a guesswork loop; always read the local B-factor and RSRZ for the residues you actually use.
- Report R-free, not R-work. R-work rewards overfitting; the R-free-minus-R-work gap is the overfitting flag, and a suspiciously small gap flags test-set leakage.
- B-factors are a within-structure relative signal only. They conflate thermal motion, static disorder, and model error, so absolute B is not portable across structures.
- pLDDT rides in the B-factor column but has opposite polarity: high value means high confidence, not high motion. Never color a predicted model "by B-factor" to infer flexibility.
- A confident predicted model can be confidently wrong. pLDDT and PAE bound the model's self-consistency, not its biological correctness; ask what context AlphaFold could not see (partner, ligand, PTM, alternative state).
- A long low-pLDDT stretch usually marks a real intrinsically disordered region, not an error. Trim it for MR/geometry pipelines, but do not conclude the protein "has no structure there" biologically.
- Cryo-EM resolution is global; local resolution varies enormously across one map. Consult the EMDB local-resolution map for peripheral domains, and keep the 0.143 (half-map) and 0.5 (map-model) FSC thresholds distinct.
- An NMR deposit is an ensemble. Validate each model and report the spread; averaging coordinates produces a physically impossible structure with distorted bonds and clashes.
- The coarse Python Ramachandran screen only decides whether to run MolProbity. Report MolProbity's numbers, which use the rama8000 reference contours the screen cannot reproduce.
- Validation precedes geometry. A distance or angle computed on an unvalidated model is only as trustworthy as the model.

## Related Skills

- structure-io - Read resolution/R-free via MMCIF2Dict and fetch the biological assembly the validation applies to
- structure-navigation - Resolve altlocs, insertion codes, and multi-model NMR files before validating per-model
- geometric-analysis - Compute the dihedrals and DSSP secondary structure this skill validates; measure only after validation passes
- structure-modification - Trim low-pLDDT residues or edit B-factors once a predicted model is validated
- alphafold-predictions - Download the AlphaFold model plus its PAE JSON that this skill reads for confidence
- modern-structure-prediction - Reconcile a re-run prediction with pLDDT/PAE/pTM when the AFDB entry is untrustworthy
- interface-analysis - Validate the assembly before interpreting an interface that only exists in it
- database-access/uniprot-access - Map validated residues back to a UniProt reference sequence

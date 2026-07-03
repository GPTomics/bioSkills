# Geometric Analysis - Usage Guide

## Overview

This skill measures geometric properties of protein structures with Biopython Bio.PDB: interatomic distances and distance matrices, bond and dihedral angles (phi/psi/chi and the Ramachandran list), superposition and RMSD, per-residue deviation, center of mass, radius of gyration, and solvent accessible surface area (SASA). Its governing idea is that no single number means "similar": every comparison metric is a choice of question, and RMSD in particular is a property of two structures GIVEN a superposition and an atom selection, not an intrinsic property of the pair. Residue contacts, contact maps, and interface/buried-surface analysis now live in structural-biology/interface-analysis; cross-protein fold comparison lives in alignment/structural-alignment; Ramachandran-as-a-quality-gate lives in structural-biology/structure-validation.

## Prerequisites

```bash
pip install biopython numpy
```

## Quick Start

Tell your AI agent what you want to do:
- "Measure the distance between residue 50 CA and residue 100 CA"
- "Superimpose these two conformations on the CA core and report the RMSD"
- "Show me where the two structures actually differ, per residue"
- "Compute per-residue relative SASA and list the buried residues"
- "Calculate phi/psi angles for every residue"

## Example Prompts

### Distances and Angles
> "What is the distance between these two atoms, and build a CA-CA distance matrix for chain A?"

> "Calculate the phi/psi angles for all residues in this structure."

> "Compute the chi1 angles for the aromatic residues."

### Superposition and RMSD
> "Superimpose the mobile structure onto the reference and report the CA RMSD."

> "The whole-molecule RMSD looks large - show me the per-residue deviation so I can see whether it is just a loop or a real difference."

> "These are two conformational states of the same protein - fit the rigid core and report core versus mobile separately."

### Surface and Shape
> "Calculate the solvent accessible surface area and tell me the probe radius used."

> "Which residues are buried versus exposed by relative SASA?"

> "What is the radius of gyration and the center of mass of this domain?"

## What the Agent Will Do

1. Parse the structure file(s) and select the requested atom set, filtering waters and heteroatoms by hetflag.
2. For distances/angles, operate on `Atom` objects (subtraction for distance) or `Vector` objects from `get_vector()` (for `calc_angle`/`calc_dihedral`).
3. For superposition, build an equal-length ordered atom correspondence, fit with `Superimposer` (SVD/Kabsch), and report `.rms` plus, when useful, the per-residue deviation that exposes outlier domination.
4. For SASA, run `ShrakeRupley.compute` at the requested level and, for burial, normalize to the Tien et al 2013 max-ASA scale rather than reporting a bare absolute area.
5. Return measurements with units and state the choices (atom selection, superposition, probe radius) that the number depends on.

## Key Functions

| Function | Purpose |
|----------|---------|
| `atom1 - atom2` | Distance between atoms |
| `calc_angle()` | Angle between 3 atoms (Vector inputs) |
| `calc_dihedral()` | Dihedral angle from 4 atoms (Vector inputs) |
| `PPBuilder().build_peptides()` | Peptides for the Ramachandran phi/psi list |
| `Superimposer` | Rigid-body superposition + RMSD (SVD/Kabsch) |
| `Bio.PDB.qcprot` | QCP superposition, faster in tight loops |
| `ShrakeRupley` | Solvent accessible surface area |

## Tips

- RMSD is not intrinsic to a pair of structures - always state the atom selection (CA? core? all-atom?) and the superposition it was measured under.
- A large global RMSD often hides a near-identical core plus one mobile loop or a hinge - use per-residue deviation to see it, and fit the rigid core separately when the biology is a domain motion.
- `Superimposer` needs equal-length ordered atom lists; it does NOT solve correspondence. Sequence-different structures need a structure-based alignment first (alignment/structural-alignment).
- To compare or rank across proteins or lengths, do not stretch RMSD - use TM-score (>0.5 = same fold) or lDDT (superposition-free); those live with alignment/structural-alignment and the AlphaFold confidence skills.
- Use `get_vector()` for angle/dihedral calculations, not `.coord`.
- A SASA number is meaningless without its probe radius (1.4A water default); absolute SASA is not portable across tools, so prefer relative SASA (Tien 2013 max-ASA) for burial.
- `Superimposer.apply` and `atom.transform` mutate coordinates in place - copy the structure first if the originals are still needed.
- Before calling two structures "different states", check the difference against B-factors, resolution, and any NMR ensemble spread; a 1-2A core RMSD is often within noise.

## Related Skills

- structure-io - Parse and write PDB/mmCIF structure files
- structure-navigation - Walk chains, residues, atoms; handle altlocs and disordered residues
- structure-modification - Transform coordinates and edit structures in place
- structural-biology/interface-analysis - Residue contacts, contact maps, and buried-surface interface analysis
- structural-biology/structure-validation - Ramachandran and omega/cis-peptide outliers as a quality gate
- alignment/structural-alignment - Cross-protein fold comparison and correspondence (TM-align, Foldseek, DALI)

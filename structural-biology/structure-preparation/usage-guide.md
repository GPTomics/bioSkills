# Structure Preparation

## Overview

A structure downloaded from the PDB or predicted by AlphaFold is not ready to dock, simulate, or run electrostatics on: X-ray models almost never resolve hydrogens, and where every proton sits depends on protonation and tautomer states that crystallographic density alone often cannot distinguish. This skill covers making a structure simulation-ready - adding hydrogens with PDBFixer or OpenMM, assigning His HID/HIE/HIP tautomers and Asn/Gln/His flips with reduce, predicting pH-dependent pKa with PROPKA, writing a PQR for electrostatics with PDB2PQR, and filling missing side-chain atoms and short loops. The through-line is that preparation ADDS inference on top of an already-inferred model, so the pH, the protonation model, and everything built must be recorded and travel with the prepared file.

## Prerequisites

```bash
pip install pdbfixer openmm
# external CLI tools installed separately:
conda install -c conda-forge reduce propka pdb2pqr
```

Familiarity with the SMCRA model and residue-id tuples (structure-navigation) and with reading structure-quality signals (structure-validation) helps decide what is safe to build.

## Quick Start

Tell your AI agent what you want to do:
- "Add hydrogens to this receptor at pH 7.4 and fill the missing side-chain atoms"
- "Assign His tautomers and Asn/Gln flips before I measure the active-site H-bonds"
- "This buried Asp - is it protonated at pH 7? Predict the pKa"
- "Make this AlphaFold model MD-ready after trimming the low-pLDDT loops"
- "Prepare this structure for docking and write a PQR for APBS"
- "Fix the missing loop and record what you built"

## Example Prompts

### Adding hydrogens and missing atoms
> "This X-ray structure has no hydrogens and a couple of truncated lysine side chains. Add hydrogens at physiological pH and complete the heavy atoms, and tell me which residues you built."

> "There is a four-residue gap in chain A between residues 112 and 117. Model the missing loop if it is short enough, and flag it as a hypothesis, not experimental."

### Protonation, tautomers, and flips
> "Before I compute interface H-bonds, resolve the Asn/Gln/His flips and pick the His tautomers by the local hydrogen-bond network."

> "My active site has a buried glutamate and a cysteine near a zinc. Do not assume standard states - predict the pKa and set the right protonation."

### Making a receptor docking/MD-ready
> "Prepare this receptor for docking: strip the cryoprotectant and buffer ions, keep the catalytic metal, add hydrogens, and assign protonation at pH 7."

> "Build an MD-ready system in the AMBER force field. Add force-field-consistent hydrogens and override any His or acid whose state the environment shifts."

### Predicted models and electrostatics
> "Prepare this ESMFold model for simulation. Trim the low-confidence regions first, then add hydrogens - do not carry the spaghetti loops into the physics."

> "I need the electrostatic potential of this protein. Assign titration states at pH 6.5 and write a PQR for APBS with the AMBER force field."

## What the Agent Will Do

1. Assess the input: resolution, missing residues/atoms, altlocs, heterogens, and whether it is experimental or predicted (predicted models need low-confidence regions trimmed first).
2. Decide what to build: short internal gaps and truncated side chains are fillable with PDBFixer; long or terminal gaps are handed to a modeling method and flagged.
3. Handle heterogens deliberately: strip buffer/cryoprotectant components, keep catalytic metals and cofactors.
4. Assign protonation: use standard states only for freely solvated residues; run PROPKA/H++ for buried, charged-clustered, metal-adjacent, or pocket residues, at an explicitly chosen pH.
5. Add hydrogens: PDBFixer or OpenMM Modeller for the bulk, reduce for Asn/Gln/His flips and His tautomer optimization by the H-bond network.
6. Emit the right output: a prepared PDB for docking/MD, or a PQR (via PDB2PQR) for Poisson-Boltzmann electrostatics with a matching force field.
7. Record provenance: the pH, the protonation model, and every atom, side chain, and loop that was built, so downstream results are reproducible.

## Tips

- Hydrogens are added, not observed: a crystal structure without hydrogens is normal, not broken. They only appear in density at roughly sub-1.2 Angstrom resolution.
- "Add H at pH 7" is a choice, not a safe default. Any titratable residue that is buried, salt-bridged, metal-adjacent, or in a pocket can deviate by several pKa units - predict it.
- His is almost never obvious: HID, HIE, and HIP give different H-bonding and charge. Let an H-bond-network optimizer (reduce) pick it.
- Asn, Gln, and His each have a 180-degree flip that density at typical resolution cannot resolve; run reduce before trusting any amide/ring H-bond.
- Missing atoms and short loops are disorder, not deletion - the atoms exist in reality. Building them is a hypothesis; never report a built loop as experimental.
- Do not let a gap-filler build long or terminal gaps; the geometry it invents is unreliable. Use a loop/homology modeler and say so.
- Strip heterogens deliberately: a blanket removal deletes catalytic metals and cofactors along with buffer salts.
- For predicted models, trim by pLDDT before preparing; pLDDT rides in the B-factor column with opposite polarity, and pocket rotamers are unreliable even where the backbone is confident.
- Match the force field and radii set across preparation and any downstream electrostatics (PDB2PQR `--ff` and the APBS run) or the charges are inconsistent.
- The prepared structure is a new object: ship it with its assumptions, or downstream H-bond, salt-bridge, interface, and docking results cannot be reproduced.

## Related Skills

- structure-validation - check resolution, altlocs, and local fit before adding inference on top
- structure-navigation - find missing residues, disorder, and altlocs; the residue-id tuple
- structure-modification - strip solvent by HETFLAG and resolve altlocs before preparation
- interface-analysis - add hydrogens here first so interface H-bond and salt-bridge geometry is meaningful
- alphafold-predictions - read pLDDT/PAE and trim low-confidence regions before preparing a predicted model
- structure-io - download the biological assembly and convert PDB/mmCIF before preparation
- chemoinformatics/virtual-screening - dock into the prepared, protonated receptor

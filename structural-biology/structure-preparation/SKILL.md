---
name: bio-structural-biology-structure-preparation
description: Prepares a deposited or predicted structure for docking, molecular dynamics, or electrostatics by adding hydrogens, assigning protonation and tautomer states, and filling missing atoms and short loops with PDBFixer, reduce, PROPKA, and PDB2PQR. Use when adding hydrogens an X-ray model never resolved; assigning His HID/HIE/HIP tautomers, Asn/Gln/His 180-degree flips, and Cys/Lys/Asp/Glu pKa-shifted protonation at a stated pH and microenvironment rather than trusting standard pKa 7; filling missing side-chain atoms and modeling short missing loops as disorder hypotheses; making a receptor docking- or MD-ready and recording what was built; preparing a predicted model after trimming low-pLDDT regions; and writing a PQR for Poisson-Boltzmann electrostatics. Keywords PDBFixer, reduce, PROPKA, PDB2PQR, protonation, tautomer, missing atoms, hydrogens, pKa, docking prep, MD prep.
tool_type: python
primary_tool: PDBFixer
---

## Version Compatibility

Reference examples tested with: pdbfixer 1.9+, openmm 8.1+

reduce/Reduce2, propka3, and pdb2pqr30 are external command-line tools installed separately (`conda install -c conda-forge reduce propka pdb2pqr` or their own pip/binary packages); the code calls them via subprocess.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Structure Preparation

**"Make this receptor simulation-ready - add hydrogens and fix the protonation"** -> add the atoms the experiment never resolved, assign pH- and environment-dependent protonation/tautomer states, and fill modeled-absent atoms and short loops.
- Python: `pdbfixer.PDBFixer` for missing atoms/residues + hydrogens at a stated pH; `openmm.app.Modeller.addHydrogens` for force-field-ready H with explicit variants
- CLI: `reduce`/`mmtbx.reduce2` for H-bond-network flip and His tautomer optimization; `propka3` for pKa prediction; `pdb2pqr30 --with-ph --ff=AMBER --titration-state-method propka` for pKa + PQR

## Governing Principle

A deposited or predicted structure is NOT simulation-ready, and "preparing" it stacks a second layer of inference on top of a model that was already inferred from data. Every atom added, every proton placed, and every loop built is a hypothesis, not a measurement - so the prepared file is a distinct object from the deposited one, and the assumptions (pH, protonation model, what was built) must travel with it or downstream H-bond, salt-bridge, interface, and docking results become unreproducible.

Hydrogens are the first trap. X-ray crystallography rarely resolves them: a hydrogen carries ~1 electron and only becomes visible in density at roughly sub-1.2 Angstrom resolution, so almost every crystal structure arrives with zero (or only polar) hydrogens. They must be ADDED, and where a hydrogen goes is decided by the PROTONATION and TAUTOMER state of its residue, which crystallographic density alone frequently cannot distinguish. His has three relevant states (HID/HIE/HIP: proton on ND1, on NE2, or both/charged); Asn, Gln, and His side-chain amide/ring groups each have a 180-degree FLIP that swaps look-alike atoms (O for N, N for C) that density at typical resolution cannot tell apart (Word et al. 1999 *J Mol Biol* 285:1735). Placing a proton "at standard pKa 7" is wrong for any residue whose microenvironment shifts its pKa: a buried Asp/Glu can stay protonated well above pH 7, a Cys in a catalytic or metal site can be a deprotonated thiolate, a Lys buried near acidic residues can lose its charge. Use a pKa predictor (PROPKA, H++) for titratable residues and an all-atom H-bond-network optimizer (reduce) for flips and His tautomers - never trust the standard state for anything not freely solvated.

Missing atoms and short missing LOOPS are the second trap. A residue truncated to Cbeta, or a chain that jumps 45 -> 58, almost always means the region was DISORDERED (too mobile to model into density), not deleted - the atoms exist in reality (see structure-navigation, structure-validation). Building them back is legitimate for making a system topologically complete, but a built side chain or loop is a GUESS among many possible conformations, must be flagged as such, and must never be reported as experimental. Long gaps, terminal extensions, and anything spanning a domain are beyond what a preparation tool should invent - hand those to a modeling method, not PDBFixer.

## Decision: tool by preparation task

| Task | Tool | Best when | Fails / misleads when |
|---|---|---|---|
| Fill missing heavy atoms, short internal loops, terminals | PDBFixer (`findMissingResidues`/`findMissingAtoms`/`addMissingAtoms`) | Truncated side chains, 1-few-residue gaps flanked by modeled residues | Long/terminal gaps, domain-scale missing regions - it builds implausible geometry |
| Add hydrogens at a pH, replace nonstandard residues (MSE, PTR) | PDBFixer `addMissingHydrogens(pH)` / `replaceNonstandardResidues` | Fast, one-call prep; standard residues in bulk-like environments | Ignores microenvironment pKa shifts and does not optimize flips/tautomers |
| Optimize Asn/Gln/His flips and His tautomer/protonation | reduce / Reduce2 (`-build` / `-FLIP`) | Resolving amide/ring orientation ambiguity by all-atom H-bond network | Treated as a pKa predictor - it optimizes geometry, not titration equilibria |
| Predict residue pKa / pH-dependent protonation | PROPKA (`propka3`), H++ | Deciding which titratable residues deviate from standard states | Reported as exact experimental pKa; empirical model, not measurement |
| Assign states + write PQR (charges/radii) for electrostatics | PDB2PQR (`pdb2pqr30 --with-ph --ff --titration-state-method propka`) | Setting up APBS/Poisson-Boltzmann; consistent charge+radius assignment | Used as a general H-adder for MD without matching the target force field |
| Add force-field-consistent H with explicit protonation variants | OpenMM `Modeller.addHydrogens(forcefield, pH, variants)` | Building an MD-ready system in a specific force field | Variants left default when a residue needs a non-standard state |
| Full MD system (solvate, neutralize, box) | OpenMM / MD prep (pointer) | After protonation is settled | Run before protonation/flips are correct - re-solvating is expensive |

## Decision: standard protonation state vs pKa predictor

| Residue / situation | Standard state at pH 7 usually fine | Use a pKa predictor + H-bond optimizer |
|---|---|---|
| Surface Asp/Glu, freely solvated | Deprotonated (-1) | Only if near a metal or H-bond partner |
| Buried or salt-bridged Asp/Glu | -- | pKa can rise several units -> may be neutral/protonated |
| Lys/Arg on the surface | Protonated (+1) | Buried Lys near acidic residues can be neutral |
| His anywhere | Ambiguous by default | Almost always: pick HID vs HIE vs HIP by local H-bonds and metal coordination |
| Cys, free | Neutral thiol | Catalytic/metal-coordinating Cys is often thiolate |
| Cys in a disulfide | No H on S (CYX) | Detect the SS bond first; do not protonate |
| Any active-site or interface residue | -- | Microenvironment dominates - predict, do not assume |

The one-line rule: standard states are defensible only for residues in a bulk-solvent-like environment; any titratable residue that is buried, charged-clustered, metal-adjacent, or in a pocket needs a predictor (PROPKA/H++) and an H-bond-network pass (reduce), with the chosen pH stated.

## PDBFixer preparation pipeline

**Goal:** Turn a raw PDB/mmCIF into a hydrogen-complete, gap-filled structure at a chosen pH, recording what was added.

**Approach:** Call the PDBFixer finders in their required order (missing residues, then nonstandard, then missing atoms) so `addMissingAtoms` sees both sets; strip crystallization heterogens while keeping (or dropping) water deliberately; then add hydrogens at an explicitly chosen pH. The order matters - hydrogens are added last, after heavy atoms exist.

```python
from pdbfixer import PDBFixer
from openmm.app import PDBFile

fixer = PDBFixer(filename='receptor.pdb')  # or PDBFixer(pdbid='1VII') to fetch from RCSB

fixer.findMissingResidues()          # short internal gaps + terminals, from SEQRES vs modeled
fixer.findNonstandardResidues()      # e.g. MSE (selenomethionine), modified residues
fixer.replaceNonstandardResidues()   # map them back to standard parents
fixer.removeHeterogens(keepWater=False)  # drop buffer ions/cryoprotectants; keepWater=True to retain
fixer.findMissingAtoms()             # truncated side chains + the residues found above
fixer.addMissingAtoms()              # build heavy atoms; built loops are HYPOTHESES, log them

# pH 7.0 is a CHOICE, not a safe default - state it and match the experimental/biological condition.
# addMissingHydrogens detects existing disulfides and leaves those Cys as CYX (no SG hydrogen).
fixer.addMissingHydrogens(pH=7.0)

with open('receptor_prepared.pdb', 'w') as out:
    PDBFile.writeFile(fixer.topology, fixer.positions, out, keepIds=True)

# missingResidues is a dict {(chain_index, residue_index): [resname, ...]} - flatten it for provenance.
built = [(ci, pos, name) for (ci, pos), names in fixer.missingResidues.items() for name in names]
print(f'built {len(fixer.missingResidues)} missing-residue segment(s), {len(built)} residue(s); pH=7.0')
```

## Optimize flips and His tautomers with reduce

**Goal:** Resolve Asn/Gln/His amide and ring orientations and His protonation by all-atom H-bond-network scoring, which density at typical resolution cannot settle.

**Approach:** Run reduce with building enabled so it adds hydrogens AND evaluates the 180-degree flip of each Asn/Gln/His plus His NH placement, choosing the orientation that optimizes the local hydrogen-bond network and minimizes clashes. reduce optimizes GEOMETRY, not titration - pair it with a pKa predictor for charge states.

```python
import subprocess

# -build runs -OH -ROTEXOH -HIS -FLIP: adds H and optimizes OH/His rotation plus Asn/Gln/His flips
# (a superset of -FLIP, not an alias). reduce scores orientation by small-probe all-atom contacts (Word 1999).
with open('receptor_reduced.pdb', 'w') as out:
    subprocess.run(['reduce', '-build', 'receptor_prepared.pdb'], stdout=out, check=True)

# Reduce2 (CCTBX/Phenix) is the maintained successor: mmtbx.reduce2 receptor_prepared.pdb
```

## Predict pKa and write a PQR for electrostatics

**Goal:** Decide which titratable residues deviate from standard states at a target pH, and emit charges+radii for Poisson-Boltzmann electrostatics.

**Approach:** PROPKA predicts per-residue pKa from the 3D environment; PDB2PQR wraps PROPKA to assign protonation at `--with-ph`, add hydrogens for the chosen force field, and write a PQR. Feed the PQR to APBS for the electrostatic potential (a separate downstream step; see the electrostatics note below).

```python
import subprocess

# propka3 writes receptor_prepared.pka; the SUMMARY lists predicted pKa vs model (standard) pKa.
subprocess.run(['propka3', 'receptor_prepared.pdb'], check=True)

# PDB2PQR assigns states at pH via PROPKA and writes charges/radii for the named force field.
# --ff must MATCH the downstream force field; --with-ph 7.0 is the stated titration condition.
subprocess.run([
    'pdb2pqr30', '--ff=AMBER', '--with-ph', '7.0',
    '--titration-state-method', 'propka', '--keep-chain',
    'receptor_prepared.pdb', 'receptor.pqr',
], check=True)
```

Electrostatics note: the PQR is the input to APBS (Jurrus et al. 2018 *Protein Sci* 27:112) for the Poisson-Boltzmann potential/surface; PDB2PQR can emit an APBS input file. Keep force field, pH, and radii set identical between preparation and the APBS run.

## Force-field-ready hydrogens with OpenMM Modeller

**Goal:** Add hydrogens consistent with a specific MD force field, forcing non-standard protonation where the environment demands it.

**Approach:** `Modeller.addHydrogens` picks the most common state per residue at the given pH and detects disulfides for Cys, but it does NOT know microenvironment pKa shifts - override with an explicit `variants` list (ASH/GLH for protonated acids, LYN for neutral Lys, HID/HIE/HIP for His) derived from a PROPKA/reduce pass. Solvation and box setup follow, at a pointer level.

```python
from openmm.app import PDBFile, Modeller, ForceField

pdb = PDBFile('receptor_prepared.pdb')
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
modeller = Modeller(pdb.topology, pdb.positions)

# variants: None per residue = let OpenMM pick at pH; override where PROPKA/reduce said otherwise.
# Set the entry for a given His to 'HID'/'HIE'/'HIP', an acid to 'ASH'/'GLH', a buried Lys to 'LYN'.
variants = modeller.addHydrogens(forcefield, pH=7.0)  # returns the chosen variant per residue

# Downstream (pointer, not this skill): modeller.addSolvent(forcefield, model='tip3p', padding=1.0*nanometer)
with open('receptor_ff_ready.pdb', 'w') as out:
    PDBFile.writeFile(modeller.topology, modeller.positions, out)
```

## Preparing a predicted (AlphaFold/ESMFold) model

**Goal:** Make a predicted model docking/MD-ready without carrying its low-confidence regions into the physics.

**Approach:** A predicted model has NO experimental hydrogens and its low-pLDDT stretches are unreliable guesses (often intrinsically disordered), so TRIM low-confidence regions FIRST, then add hydrogens/protonation. pLDDT rides in the B-factor column (opposite polarity to a real B-factor); use it to cut, not to color as mobility (see alphafold-predictions). Pocket rotamers are the least reliable atoms even where backbone pLDDT is high, so verify the binding site before docking.

```python
from pdbfixer import PDBFixer
from openmm.app import PDBFile

# 1) Trim low-pLDDT residues (pLDDT<50-70 = unreliable) BEFORE preparation; phenix.process_predicted_model
#    does this + a PAE domain split for MR. Here: a minimal B-factor(=pLDDT) filter as illustration.
# 2) Then run the PDBFixer pipeline above on the trimmed model to add H and any missing side-chain atoms.
fixer = PDBFixer(filename='af_model_trimmed.pdb')
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(pH=7.0)  # predicted models never carry experimental hydrogens
with open('af_model_prepared.pdb', 'w') as out:
    PDBFile.writeFile(fixer.topology, fixer.positions, out, keepIds=True)
```

## Common Errors

| Symptom | Cause | Fix |
|---|---|---|
| Prepared file has no hydrogens | X-ray models rarely resolve H; parsing does not add them | Run `addMissingHydrogens(pH=...)` or `Modeller.addHydrogens` explicitly |
| His H-bonds/metal coordination look wrong | Default HIE/HID guessed without the local network; wrong tautomer | Let reduce pick the neutral tautomer (HID vs HIE) by H-bond network; use PROPKA only to decide the CHARGE state (HIP vs neutral) - they answer different questions, so reconcile rather than pick one |
| Buried Asp/Glu deprotonated but should be neutral | Standard pKa 7 assumed; buried pKa is shifted up | Predict pKa (PROPKA/H++); protonate the residue (ASH/GLH) |
| Catalytic Cys modeled as neutral thiol | Standard state assumed in a metal/active site | Predict pKa / check metal coordination; set thiolate or CYX for disulfides |
| Asn/Gln side chain H-bonds backwards | 180-degree amide flip not resolved (O/N indistinguishable in density) | Run reduce with flips enabled before analysis |
| `addMissingAtoms` builds a wild loop | A long/terminal gap handed to a gap-filler that only does short loops | Do not build long gaps here; use a loop/homology modeler and flag it |
| Catalytic metal or cofactor gone after prep | `removeHeterogens()` stripped all non-water heterogens | Keep needed heterogens: filter deliberately, do not blanket-remove |
| Missing residues not built | `findMissingAtoms` called before `findMissingResidues` | Call finders in order: residues, nonstandard, then atoms |
| pKa/protonation differs from a paper | Different pH or predictor; states are pH- and method-dependent | State the pH and tool; treat predicted pKa as a model, not a measurement |
| Predicted pKa sits close to the working pH | The protonation state is genuinely ambiguous, and empirical predictors are weakest at metal and strongly-coupled active sites | Test both states (or run constant-pH MD); at metal/catalytic sites treat the predicted pKa as a weak prior and cross-check coordination geometry/literature |
| Downstream results not reproducible | Prepared file shipped without its assumptions | Record pH, protonation model, and every built atom/loop as provenance |
| APBS charges look wrong | PDB2PQR `--ff` did not match the downstream force field/radii | Match `--ff` and radii set across preparation and APBS |
| Predicted model docks into a garbage pocket | Low-pLDDT/rotamer-unreliable region kept, or wrong apo/holo state | Trim by pLDDT first; verify the pocket conformation before docking |

## Related Skills

- structure-validation - check resolution, altlocs, and the region of interest before adding inference on top
- structure-navigation - identify missing residues, disorder, altlocs, and the (hetflag, resseq, icode) id tuple
- structure-modification - strip solvent by HETFLAG and resolve altlocs before preparation; never overwrite pLDDT-in-B
- interface-analysis - add hydrogens here first so H-bond and salt-bridge geometry across an interface is meaningful
- alphafold-predictions - read pLDDT/PAE and trim low-confidence regions before preparing a predicted model
- structure-io - download the biological assembly and convert PDB/mmCIF before preparation
- chemoinformatics/virtual-screening - dock into the prepared, protonated receptor

## References

- Eastman P, Swails J, Chodera JD, et al. 2017. OpenMM 7: rapid development of high performance algorithms for molecular dynamics. *PLoS Comput Biol* 13(7):e1005659. doi:10.1371/journal.pcbi.1005659 (PDBFixer ships with OpenMM)
- Word JM, Lovell SC, Richardson JS, Richardson DC. 1999. Asparagine and glutamine: using hydrogen atom contacts in the choice of side-chain amide orientation. *J Mol Biol* 285(4):1735-1747. doi:10.1006/jmbi.1998.2401 (reduce Asn/Gln/His flips)
- Olsson MHM, Sondergaard CR, Rostkowski M, Jensen JH. 2011. PROPKA3: consistent treatment of internal and surface residues in empirical pKa predictions. *J Chem Theory Comput* 7(2):525-537. doi:10.1021/ct100578z
- Sondergaard CR, Olsson MHM, Rostkowski M, Jensen JH. 2011. Improved treatment of ligands and coupling effects in empirical calculation and rationalization of pKa values. *J Chem Theory Comput* 7(7):2284-2295. doi:10.1021/ct200133y
- Dolinsky TJ, Nielsen JE, McCammon JA, Baker NA. 2004. PDB2PQR: an automated pipeline for the setup of Poisson-Boltzmann electrostatics calculations. *Nucleic Acids Res* 32:W665-W667. doi:10.1093/nar/gkh381
- Jurrus E, Engel D, Star K, et al. 2018. Improvements to the APBS biomolecular solvation software suite. *Protein Sci* 27(1):112-128. doi:10.1002/pro.3280
- Anandakrishnan R, Aguilar B, Onufriev AV. 2012. H++ 3.0: automating pK prediction and the preparation of biomolecular structures for atomistic molecular modeling and simulations. *Nucleic Acids Res* 40:W537-W541. doi:10.1093/nar/gks375

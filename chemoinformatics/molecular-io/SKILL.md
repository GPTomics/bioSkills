---
name: bio-molecular-io
description: Reads, writes, and converts molecular file formats (SMILES, InChI, SDF V2000/V3000, MOL2, PDB, MMTF) using RDKit and Open Babel with rigorous handling of aromaticity perception, stereochemistry, implicit/explicit hydrogens, kekulization, and salt/fragment separation. Use when loading chemical libraries, debugging parse failures, or preparing molecules for downstream standardization, descriptor calculation, or docking.
tool_type: python
primary_tool: RDKit
---

## Version Compatibility

Reference examples tested with: RDKit 2024.09+, Open Babel 3.1.1+, ChEMBL structure_pipeline 1.2+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `obabel -V`; `obabel -L formats`

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Molecular I/O

Parse, write, and convert molecular file formats. Most downstream errors trace back to silent I/O issues: incorrect aromaticity perception, lost stereochemistry, mishandled charges, dropped stereo bonds, or non-canonical tautomers. This skill enumerates each format's failure modes and prescribes the correct toolchain for each scenario.

For full standardization (canonicalization, salt stripping, tautomer enumeration) see `chemoinformatics/molecular-standardization`. For generating 3D conformers from parsed 2D molecules, see `chemoinformatics/conformer-generation`.

## Format Taxonomy

| Format | Dim | Stereo | Charges | Strength | Fails when |
|--------|-----|--------|---------|----------|------------|
| SMILES | 2D | Explicit `/\@` | Net charges only | Compact, web-friendly, fast parse | Loses absolute coordinates; aromatic perception ambiguous across toolkits; tautomers not canonical |
| InChI | 2D | Layer 't,s' | pH-fixed std layer 'p' | Canonical by construction; cross-toolkit identity | Loses tautomers (std InChI merges); loses stereo at C^x metals; large molecules timeout |
| SDF V2000 | 2D/3D | Wedge bonds | M CHG line | Industry default; metadata via tags | 999-atom limit; cannot encode multi-component reactions; query atoms ambiguous |
| SDF V3000 | 2D/3D | Wedge + stereo flag | Inline charge | No atom limit; query support; rich properties | Some software (legacy) cannot read; verbose |
| MOL2 (Tripos) | 3D | Wedge bonds | Per-atom partial | SYBYL atom types preserved for docking | Atom-type dialects diverge (SYBYL vs Corina); RDKit MOL2 parser brittle |
| PDB | 3D | None | None standard | Universal protein format | No bond orders; aromatic perception lost; ligand names truncated to 3 chars |
| PDBQT | 3D | None | Gasteiger / AD4 | AutoDock-ready; torsion tree encoded | Specific to docking; no aromaticity layer |
| BinaryCIF (MMTF retired) | 3D | Encoded | Encoded | BinaryCIF (`.bcif`) is the current compact structural format; MMTF was retired by RCSB in July 2024 (read-only legacy) | Not all toolkits parse; binary format |
| CDX/CDXML | 2D | Drawing | Drawing | ChemDraw native | Not a structural format; converts unreliably |
| InChIKey | Hash | Stereo layer | n/a | Database key, fast lookup | Hash collisions ~10^-9 but possible; cannot recover structure |

## Aromaticity Perception (most common silent error)

Different toolkits perceive aromaticity differently. The same SMILES round-tripped between toolkits may produce different canonical strings and different fingerprints.

| Model | Toolkit | Rule | Symptom of mismatch |
|-------|---------|------|---------------------|
| Daylight | OpenEye, Daylight | 4n+2 π on planar ring | Furan, thiophene aromatic |
| RDKit default | RDKit | Daylight-like with extensions for fused / N-containing | Compatible with Daylight for drug-like molecules |
| MDL | Indigo, ChemAxon | Reduced (only pyrrole-type) | Pyrrole NH aromatic but tropone non-aromatic |
| OpenEye | OEAroModel | Several modes | Charged thiophene non-aromatic in MDL but aromatic in OpenEye |

**Fix:** Always re-canonicalize via the toolkit doing analysis. Never trust a SMILES produced by toolkit A as input to toolkit B without `SetAromaticity(rdkit.Chem.AromaticityModel)`.

## Stereochemistry Layers

Stereo loss is the second most common silent error. Each format encodes stereo differently:

- SMILES: `@/@@` for tetrahedral, `/` and `\` for cis/trans double bonds
- InChI: separate `/t` (tetrahedral), `/s` (stereo flag), `/m` (mirror) layers
- SDF: wedge/hash bond + parity 0/1/2; cis/trans encoded via bond direction
- MOL2: explicit stereo flag

**Round-trip tests:** If `Chem.MolToSmiles(Chem.MolFromSmiles(smi))` does not preserve `@` and `/\`, the molecule was sanitized without `Chem.RemoveStereochemistry()`. If `MolFromMolFile` returns a molecule missing wedge bonds, the SDF used parity-only encoding (legacy).

## Reading SMILES with Stereo Preservation

**Goal:** Parse SMILES while preserving stereo and aromatic-flag consistency.

**Approach:** Use `Chem.MolFromSmiles(smi)` with sanitization on, verify with round-trip canonicalization, and set explicit stereochemistry where the toolkit's perception missed it.

```python
from rdkit import Chem
from rdkit.Chem import AllChem

def parse_smiles_safe(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None, 'parse_failure'
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    canon = Chem.MolToSmiles(mol)
    round_trip = Chem.MolFromSmiles(canon)
    if Chem.MolToSmiles(round_trip) != canon:
        return mol, 'round_trip_unstable'
    return mol, 'ok'
```

## Reading SDF with Property Carryover

**Goal:** Load a multi-record SDF preserving per-molecule properties (Name, ID, IC50, etc.) used by downstream filtering and ML labeling.

**Approach:** Iterate via `SDMolSupplier(removeHs=False, sanitize=True)`, filter `None` (parse failures), and capture properties via `mol.GetPropsAsDict()`.

```python
from rdkit import Chem

supplier = Chem.SDMolSupplier('library.sdf', removeHs=False, sanitize=True)
mols = []
fails = []
for i, mol in enumerate(supplier):
    if mol is None:
        fails.append(i)
        continue
    props = mol.GetPropsAsDict()
    mols.append((mol, props))
print(f'parsed: {len(mols)}; failed: {len(fails)}')
```

If a large fraction fails, try `sanitize=False` then `Chem.SanitizeMol(mol, catchErrors=True)` to identify per-step failures (kekulization, valence, aromaticity).

## Open Babel for MOL2 / PDBQT

RDKit's MOL2 parser is incomplete (SYBYL atom-type sets differ). Open Babel is more robust for MOL2 and PDBQT.

```python
from openbabel import pybel

mols = list(pybel.readfile('mol2', 'ligands.mol2'))
for mol in mols:
    smi = mol.write('smi').strip().split()[0]
    inchi = mol.write('inchi').strip()
```

For docking output PDBQT, use Open Babel rather than RDKit:
```python
import subprocess
subprocess.run(['obabel', 'docked.pdbqt', '-O', 'docked.sdf'], check=True)
```

## InChI for Canonical Identity

InChI is the only canonical 2D representation guaranteed across toolkits. Standard InChI ignores tautomers and metal stereo; non-standard variants preserve them with flags.

```python
from rdkit.Chem.inchi import MolToInchi, MolToInchiKey, InchiToInchiKey

mol = Chem.MolFromSmiles('c1ccc2c(c1)cccc2')
inchi = MolToInchi(mol)
key = MolToInchiKey(mol)

inchi_fixedH, retcode, msg, _, _ = Chem.MolToInchiAndAuxInfo(mol, options='/FixedH')
```

**Caveat:** Two molecules with identical std InChI may be different tautomers. Use `/FixedH` for tautomer-distinguishing InChI when needed.

## Per-Format Failure Modes

### SMILES -- ambiguous aromaticity

**Trigger:** Input from non-RDKit source (ChemAxon, OpenEye, Daylight) round-tripping into RDKit.

**Mechanism:** RDKit perceives aromaticity on input. Aromatic flags from origin toolkit are overwritten.

**Symptom:** Fingerprints differ between toolkits for "identical" molecules; database joins by canonical SMILES miss records.

**Fix:** Always re-canonicalize within the analysis toolkit. For cross-toolkit identity, use InChIKey not canonical SMILES.

### SDF V2000 -- atom count >999

**Trigger:** Large molecules (peptides, oligonucleotides, dendrimers).

**Mechanism:** V2000 header uses fixed 3-character atom count field.

**Symptom:** Truncated atom block; parse failure with cryptic error.

**Fix:** Switch to V3000: `Chem.SDWriter('out.sdf', forceV3000=True)`. RDKit auto-detects V3000 on read; explicitly request on write.

### SDF -- wedge bond orientation lost

**Trigger:** SDF written by tools that use parity flags only (older ISIS-Draw, some pipeline tools).

**Mechanism:** Parity alone is ambiguous without geometric coordinates; RDKit reads parity but cannot re-render wedges.

**Symptom:** Drawn molecule shows undefined stereo despite SDF carrying parity bits.

**Fix:** After read, `Chem.AssignStereochemistryFrom3D(mol)` if 3D coords present; otherwise stereo must be re-derived from SMILES with wedges.

### PDB ligand -- no bond orders

**Trigger:** Parsing ligand from PDB entry (e.g., extracting co-crystal ligand).

**Mechanism:** PDB stores only atoms + CONECT; bond orders inferred by RDKit's `AssignBondOrdersFromTemplate` which requires a template molecule.

**Symptom:** All bonds single; aromatic rings non-aromatic; valences wrong.

**Fix:** Use `AllChem.AssignBondOrdersFromTemplate(template, ligand)` where `template` is a SMILES-derived mol of the expected ligand structure. Or use the PDB Ligand Expo SDF.

### MOL2 -- SYBYL atom type dialect

**Trigger:** MOL2 produced by Corina, MOE, or Schrodinger.

**Mechanism:** SYBYL atom types are not perfectly standardized across vendors; RDKit's parser handles canonical SYBYL.

**Symptom:** Mol returns as `None` or with wrong atom types (`Cl` vs `Cl.O` peroxide-style).

**Fix:** Convert via Open Babel as intermediate: `obabel input.mol2 -O temp.sdf` then read SDF.

### Open Babel pybel -- import path

**Trigger:** Code written for Open Babel 2.x.

**Mechanism:** OB 3.x reorganized: `import pybel` no longer works.

**Symptom:** `ModuleNotFoundError: No module named 'pybel'`.

**Fix:** `from openbabel import pybel`.

## Charge Models on I/O

| Source | Charges in file | Use for |
|--------|-----------------|---------|
| Parsed SMILES | Net formal charges only | Storage, similarity, ML training |
| Parsed PDB | Atomic charges typically absent | Always re-assign for downstream |
| `obabel --partialcharge gasteiger` | Gasteiger Marsili (empirical) | AutoDock Vina, fast |
| AM1-BCC (AmberTools antechamber) | Semi-empirical | MD, FEP setup |
| RESP (psi4, Gaussian) | Quantum ESP-fitted | High-accuracy MD, FEP |

The charge model **must** match the downstream method. Mixing AM1-BCC ligand charges with TIP3P water + AMBER protein is valid; Gasteiger charges are unsuitable for MD.

## Drawing for QC

Always draw a random subset of parsed molecules. Wrong stereo, missing rings, and broken aromaticity show immediately.

```python
from rdkit.Chem.Draw import rdMolDraw2D

def draw_grid(mols, fname, mols_per_row=5, sub_img_size=(250, 200)):
    from rdkit.Chem.Draw import MolsToGridImage
    img = MolsToGridImage(mols[:25], molsPerRow=mols_per_row, subImgSize=sub_img_size,
                          legends=[m.GetProp('_Name') if m.HasProp('_Name') else ''
                                   for m in mols[:25]])
    img.save(fname)
```

`MolsToGridImage` returns PIL image; for headless servers use `MolDraw2DCairo` directly.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| `Chem.MolFromSmiles` returns None | Invalid SMILES, bad parentheses, ring not closed | Try `sanitize=False`, inspect with `Chem.MolFromSmiles(smi, sanitize=False)` |
| Round-trip SMILES changes | Aromaticity perception drift | Always canonicalize within analysis toolkit |
| All bonds single in PDB ligand | PDB has no bond orders | `AllChem.AssignBondOrdersFromTemplate(template, mol)` |
| Stereo lost on SDF write | `Chem.RemoveStereochemistry` called | Use `MolToMolBlock(mol, kekulize=True)` |
| MOL2 parse returns None | RDKit MOL2 parser incomplete for vendor dialects | Convert via Open Babel intermediate |
| InChI differs for "same" molecule | Different tautomers, charges, or stereo | Use `/FixedH` to retain tautomer; compare without standardization |
| Fingerprints differ across toolkits | Aromaticity model difference | Use InChIKey for identity; re-canonicalize for similarity |

## References

- Heller et al., *J. Cheminformatics* 7:23 (2015) -- InChI v1.05 specification.
- O'Boyle, *J. Cheminformatics* 4:22 (2012) -- canonical SMILES algorithms across toolkits.
- Bento et al., *J. Cheminformatics* 12:51 (2020) -- ChEMBL structure pipeline.
- Riniker & Landrum, *J. Chem. Inf. Model.* 55:2562-2574 (2015) -- ETKDG conformer embedding (RDKit knowledge-distance-geometry).

## Related Skills

- chemoinformatics/molecular-standardization - Salt stripping, tautomer canonicalization, neutralization
- chemoinformatics/molecular-descriptors - Calculate fingerprints and properties from parsed molecules
- chemoinformatics/conformer-generation - Generate 3D coordinates from 2D inputs
- chemoinformatics/virtual-screening - Prepare ligands for docking
- structural-biology/structure-io - Protein structure handling (PDB, mmCIF)

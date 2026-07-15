---
name: bio-molecular-io
description: Reads, writes, and converts molecular file formats (SMILES, InChI, SDF V2000/V3000, MOL2, PDB, and BinaryCIF) using RDKit and Open Babel with rigorous handling of aromaticity perception, stereochemistry, implicit/explicit hydrogens, kekulization, and salt/fragment separation. Use when loading chemical libraries, debugging parse failures, or preparing molecules for downstream standardization, descriptor calculation, or docking.
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
| SMILES | 2D | Atom chirality `@/@@`; double-bond `/` and `\` | Atom-local formal charges | Compact, web-friendly, fast parse | Loses absolute coordinates; aromatic perception ambiguous across toolkits; tautomers not canonical |
| InChI | 2D | `/b`, `/t`, `/m`, `/s` stereo sublayers | `/q` charge and `/p` added/removed-proton sublayers; `/p` is not a pH model | Canonical by construction; cross-toolkit identity | Standard InChI normalizes mobile-H forms; limited organometallic stereo; large molecules may require special handling |
| SDF V2000 | 2D/3D | Wedge bonds | M CHG line | Industry default; metadata via tags | 999-atom limit; cannot encode multi-component reactions; query atoms ambiguous |
| SDF V3000 | 2D/3D | Wedge + stereo flag | Inline charge | No atom limit; query support; rich properties | Some software (legacy) cannot read; verbose |
| MOL2 (Tripos) | 3D | Common records rely on 3D coordinates and toolkit perception; no portable explicit stereo field | Per-atom partial | SYBYL atom types preserved for docking | Atom-type dialects diverge (SYBYL vs Corina); RDKit MOL2 parser brittle |
| PDB | 3D | None | None standard | Universal protein format | No bond orders; aromatic perception lost; ligand names truncated to 3 chars |
| PDBQT | 3D | None | Gasteiger / AD4 | AutoDock-ready; torsion tree encoded | Specific to docking; no aromaticity layer |
| BinaryCIF (MMTF retired) | 3D | Encoded | Encoded | BinaryCIF (`.bcif`) is the current compact structural format; RCSB stopped serving MMTF files on July 2, 2024 and recommends BinaryCIF (RCSB PDB 2024) | Not all toolkits parse; binary format |
| CDX/CDXML | 2D | Drawing | Drawing | ChemDraw native | Not a structural format; converts unreliably |
| InChIKey | Hash | Stereo layer | n/a | Database key, fast lookup | Collision probability depends on key block and collection size; cannot recover structure |

## Aromaticity Perception (most common silent error)

Different toolkits perceive aromaticity differently. The same SMILES round-tripped between toolkits may produce different canonical strings and different fingerprints.

| Model | Toolkit | Rule | Symptom of mismatch |
|-------|---------|------|---------------------|
| Daylight | OpenEye, Daylight | 4n+2 π on planar ring | Furan, thiophene aromatic |
| RDKit default | RDKit | Daylight-like with extensions for fused / N-containing | Compatible with Daylight for drug-like molecules |
| MDL | Available in several toolkits, including RDKit as `AROMATICITY_MDL` | Five-membered rings are not aromatic unless part of a fused aromatic system; only C/N and one-electron donors qualify; exocyclic double bonds exclude an atom | Five-membered heteroaromatics and exocyclic-bond systems may differ from the default RDKit model |
| OpenEye | OEAroModel | Several modes | Charged thiophene non-aromatic in MDL but aromatic in OpenEye |

**Fix:** Always re-canonicalize via the toolkit doing analysis. If aromaticity must be reassigned explicitly in RDKit, use a concrete model, for example `Chem.SetAromaticity(mol, Chem.AromaticityModel.AROMATICITY_RDKIT)`, after the molecule is in an appropriate sanitized or kekulized state.

## Stereochemistry Layers

Stereo loss is the second most common silent error. Each format encodes stereo differently:

- SMILES: `@/@@` for tetrahedral, `/` and `\` for cis/trans double bonds
- InChI: `/b` for double-bond stereo, `/t` for tetrahedral stereo, and `/m` plus `/s` for inversion/overall stereo type
- SDF: wedge/hash bond + parity 0/1/2; cis/trans encoded via bond direction
- MOL2: common Tripos records provide 3D coordinates but no portable wedge or explicit stereo field; verify toolkit perception by round trip

**Round-trip tests:** If `Chem.MolToSmiles(Chem.MolFromSmiles(smi), isomericSmiles=True)` does not preserve the represented stereochemistry, inspect whether the source contained stereo markers and whether any step called `Chem.RemoveStereochemistry()` or discarded stereochemical coordinates/bond directions. Sanitization alone does not intentionally remove valid stereochemistry. If `MolFromMolFile` returns a molecule missing wedge bonds, inspect the source's coordinates, bond directions, and parity encoding.

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

InChI is a standardized, canonical structure identifier designed for cross-database and cross-toolkit interoperability (Heller et al. 2015; O'Boyle 2012). Standard InChI normalizes many mobile-hydrogen tautomers and has limitations for some metal stereochemistry; non-standard options such as `/FixedH` can distinguish additional representations. InChIKey is a fixed-length hash, so use full InChI or standardized structures when a suspected collision must be resolved (InChI Trust technical FAQ).

```python
from rdkit.Chem.inchi import MolToInchi, MolToInchiKey, InchiToInchiKey

mol = Chem.MolFromSmiles('c1ccc2c(c1)cccc2')
inchi = MolToInchi(mol)
key = MolToInchiKey(mol)

inchi_fixedH, aux_info = Chem.MolToInchiAndAuxInfo(mol, options='/FixedH')
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

**Fix:** Switch to V3000. RDKit auto-detects V3000 on read; explicitly request it on the writer:

```python
writer = Chem.SDWriter('out.sdf')
writer.SetForceV3000(True)
writer.write(mol)
writer.close()
```

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
| Parsed SMILES | Atom-local formal charges; no partial-charge model | Storage, similarity, ML training |
| Parsed PDB | Atomic charges typically absent | Always re-assign for downstream |
| `obabel --partialcharge gasteiger` | Gasteiger-Marsili partial charges (empirical) | Workflows that explicitly require Gasteiger charges; Vina/Vinardo scoring itself does not require assigned atom charges |
| AM1-BCC (AmberTools antechamber) | Semi-empirical | MD, FEP setup |
| RESP (psi4, Gaussian) | Restrained fit to a quantum-mechanical ESP; protocol-specific | Force-field workflows parameterized for that RESP protocol |

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
| Stereo lost on SDF write | Stereo was absent, removed, or not represented by coordinates/bond directions | Verify assigned chiral tags and bond stereo before writing; preserve suitable 2D/3D coordinates and inspect the round trip |
| MOL2 parse returns None | RDKit MOL2 parser incomplete for vendor dialects | Convert via Open Babel intermediate |
| InChI differs for "same" molecule | Different tautomers, charges, or stereo | Use `/FixedH` to retain tautomer; compare without standardization |
| Fingerprints differ across toolkits | Aromaticity model difference | Use InChIKey for identity; re-canonicalize for similarity |

## References

- Heller et al., *J. Cheminformatics* 7:23 (2015) -- InChI design, layout, and algorithms; software version 1.04. https://doi.org/10.1186/s13321-015-0068-4
- O'Boyle, *J. Cheminformatics* 4:22 (2012) -- Universal SMILES representation based on InChI. https://doi.org/10.1186/1758-2946-4-22
- InChI Trust, "Technical FAQ" -- InChIKey layout and collision estimates. https://www.inchi-trust.org/technical-faq/
- Daylight Chemical Information Systems, SMILES theory -- atom-local formal-charge and stereochemical syntax. https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html
- RDKit Book, "Aromaticity" -- RDKit and MDL aromaticity-model rules. https://www.rdkit.org/docs/RDKit_Book.html#aromaticity
- RCSB PDB (2024), "Removal of MMTF files from RCSB PDB." https://www.rcsb.org/news/6661c73362451e4e35915f7b

## Related Skills

- chemoinformatics/molecular-standardization - Salt stripping, tautomer canonicalization, neutralization
- chemoinformatics/molecular-descriptors - Calculate fingerprints and properties from parsed molecules
- chemoinformatics/conformer-generation - Generate 3D coordinates from 2D inputs
- chemoinformatics/virtual-screening - Prepare ligands for docking
- structural-biology/structure-io - Protein structure handling (PDB, mmCIF)

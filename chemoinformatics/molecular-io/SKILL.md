---
name: bio-molecular-io
description: Reads, writes, and converts molecular file formats (SMILES, SDF, MOL2, PDB) using RDKit and Open Babel. Handles structure parsing, canonicalization, and full standardization pipeline including sanitization, normalization, and tautomer canonicalization. Use when loading chemical libraries, converting formats, or preparing molecules for analysis.
tool_type: python
primary_tool: RDKit
---

# Molecular I/O

Read, write, and convert molecular file formats with structure standardization.

## Supported Formats

| Format | Extension | Use Case |
|--------|-----------|----------|
| SMILES | .smi | Text representation, databases |
| SDF/MOL | .sdf, .mol | 3D structures, compound libraries |
| MOL2 | .mol2 | Docking, force field atoms |
| PDB | .pdb | Protein-ligand complexes |

## Reading Molecules

```python
from rdkit import Chem
from rdkit.Chem import AllChem

# From SMILES
mol = Chem.MolFromSmiles('CCO')

# From SDF file (single molecule)
mol = Chem.MolFromMolFile('molecule.mol')

# From SDF file (multiple molecules)
supplier = Chem.SDMolSupplier('library.sdf')
molecules = [mol for mol in supplier if mol is not None]
print(f'Loaded {len(molecules)} molecules')

# From SMILES file
with open('compounds.smi') as f:
    molecules = []
    for line in f:
        parts = line.strip().split()
        if parts:
            mol = Chem.MolFromSmiles(parts[0])
            if mol:
                mol.SetProp('_Name', parts[1] if len(parts) > 1 else '')
                molecules.append(mol)
```

## Writing Molecules

```python
from rdkit import Chem

# To SMILES
smiles = Chem.MolToSmiles(mol)  # Canonical SMILES
smiles_iso = Chem.MolToSmiles(mol, isomericSmiles=True)  # With stereochemistry

# To SDF file
writer = Chem.SDWriter('output.sdf')
for mol in molecules:
    writer.write(mol)
writer.close()

# To MOL block (string)
mol_block = Chem.MolToMolBlock(mol)
```

## Structure Standardization

Use rdMolStandardize module (Python MolStandardize was removed Q1 2024).

```python
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

def standardize_molecule(mol):
    '''
    Full standardization pipeline.
    Order: Sanitize -> Normalize -> Neutralize -> Canonicalize tautomer -> Strip salts
    '''
    if mol is None:
        return None

    # Sanitize (assign valences, kekulize)
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        return None

    # Normalize (standardize functional groups)
    normalizer = rdMolStandardize.Normalizer()
    mol = normalizer.normalize(mol)

    # Neutralize charges where possible
    uncharger = rdMolStandardize.Uncharger()
    mol = uncharger.uncharge(mol)

    # Canonicalize tautomers
    enumerator = rdMolStandardize.TautomerEnumerator()
    mol = enumerator.Canonicalize(mol)

    # Remove salts/fragments (keep largest)
    remover = rdMolStandardize.FragmentRemover()
    mol = remover.remove(mol)

    return mol

# Standardize a library
standardized = [standardize_molecule(m) for m in molecules]
standardized = [m for m in standardized if m is not None]
```

## Open Babel Conversion

For format conversions not supported by RDKit.

```python
# Open Babel 3.x import (not 'import pybel')
from openbabel import pybel

# Read MOL2 (better supported in Open Babel)
mols = list(pybel.readfile('mol2', 'ligands.mol2'))

# Convert to SDF
output = pybel.Outputfile('sdf', 'output.sdf', overwrite=True)
for mol in mols:
    output.write(mol)
output.close()

# Format conversion
for mol in pybel.readfile('pdb', 'complex.pdb'):
    mol.write('mol2', 'ligand.mol2', overwrite=True)
```

## Molecular Drawing

Use rdMolDraw2D (legacy Draw.MolToImage deprecated).

```python
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

def draw_molecule(mol, filename, size=(400, 300)):
    '''Draw molecule to PNG file.'''
    drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    with open(filename, 'wb') as f:
        f.write(drawer.GetDrawingText())

# Draw with highlighting
def draw_with_substructure(mol, pattern, filename):
    '''Highlight substructure match.'''
    match = mol.GetSubstructMatch(Chem.MolFromSmarts(pattern))
    drawer = rdMolDraw2D.MolDraw2DCairo(400, 300)
    drawer.DrawMolecule(mol, highlightAtoms=match)
    drawer.FinishDrawing()
    with open(filename, 'wb') as f:
        f.write(drawer.GetDrawingText())
```

## Related Skills

- molecular-descriptors - Calculate properties after loading
- similarity-searching - Compare loaded molecules
- virtual-screening - Prepare ligands for docking

---
name: bio-structural-biology-structure-modification
description: Modifies protein structures in place with Biopython Bio.PDB - transforms coordinates, strips waters/heteroatoms, overloads the B-factor column, renumbers, and builds entities. Use when applying a rotation matrix and needing to know whether it is row-convention (Entity.transform, Superimposer) or column-convention (REMARK 350 / _pdbx_struct_oper_list assembly operators) so geometry is not silently mirrored; when overloading B-factors with pLDDT/conservation for coloring and needing to preserve the destroyed originals; when stripping solvent by HETFLAG (r.id[0]) rather than residue name so catalytic metals and cofactors survive; and when building or copying entities through StructureBuilder/Select without breaking SMCRA parent-child links or the (hetflag, resseq, icode) id tuple. Keywords transform, rotation matrix, occupancy, assembly operators.
tool_type: python
primary_tool: Bio.PDB
goal_approach_exempt: true
---

## Version Compatibility

Reference examples tested with: biopython 1.83+, numpy 1.26+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Structure Modification

**"Move this chain onto that one and strip the waters"** -> mutate coordinates and the entity tree in place, then write a new file.
- Python: `Entity.transform(rot, tran)` for coordinates, `detach_child` / `PDBIO(select=...)` for filtering, `StructureBuilder` for building

## Governing Principle: every edit mutates in place, and the rotation convention is the load-bearing trap

Bio.PDB has no immutable copy semantics. `atom.coord = ...`, `residue.id = ...`, `chain.detach_child(...)`, and `Entity.transform(...)` all mutate the parsed object directly, so the moment a downstream step still needs the original, a `copy.deepcopy` must be taken first (a plain reference is not a copy).

The trap that silently corrupts geometry is the rotation convention. Bio.PDB `Superimposer`, `SVDSuperimposer`, and `Entity.transform(rot, tran)` apply the transform as `dot(coords, rot) + tran` - coordinates are treated as ROW vectors post-multiplied by `rot`, so the `rot` these classes hand back is the TRANSPOSE of the textbook rotation matrix. Biological-assembly operators are the opposite: REMARK 350 and mmCIF `_pdbx_struct_oper_list` matrices are COLUMN-convention (`R @ x + t`). Feeding a column-convention `R` straight into `Entity.transform` (or writing `np.dot(R, atom.coord)` against a row-convention source) applies the transpose and yields a mirrored or wrongly-rotated structure that still looks plausible. Prefer `Entity.transform` / `atom.transform` (which own the row convention) over hand-rolled `np.dot`, and transpose any column-convention operator before passing it in.

Three more edits destroy data quietly: overloading the B-factor column with a per-residue scalar (pLDDT, conservation) DESTRUCTIVELY overwrites the real temperature factors - and for AlphaFold models the column already IS pLDDT, so overwrite it and the confidence signal is gone; save the originals first. Stripping solvent by residue NAME instead of the HETFLAG (`r.id[0]`) deletes functional metals, cofactors, and modified residues (MSE) mid-chain. And building or copying entities without wiring the SMCRA parent-child links, or renumbering without carrying the full `(hetflag, resseq, icode)` id tuple, makes the writer emit broken or collided records.

## Decision: which transform path

| Matrix source | Convention | Apply as | Failure if mixed |
|---|---|---|---|
| `Superimposer.rotran` / `SVDSuperimposer.get_rotran` | row (`coords @ rot`) | `Entity.transform(rot, tran)` | none - same convention |
| `Entity.transform` / `atom.transform` | row (`coords @ rot`) | pass `rot` as-is | none |
| REMARK 350 / `_pdbx_struct_oper_list` assembly operators | column (`R @ x + t`) | `Entity.transform(R.T, t)` | column `R` applied row -> mirrored/rotated wrong |
| `Bio.PDB.vectors.rotaxis(theta, Vector)` | row (built for `.transform`) | `Entity.transform(rot, tran)` | none |
| Raw math / textbook `R` via `np.dot` | column (`R @ x`) | `R @ coord + t` explicitly, consistently | inconsistent left/right multiply |

## Decision: how to strip solvent and hetero

| Strategy | Filter | Deletes | Use when |
|---|---|---|---|
| By HETFLAG, water only | `r.id[0] == 'W'` | ordered/crystallographic waters | safe default before docking/MD prep |
| By explicit deny-list | `r.resname in {'HOH','SO4','GOL','EDO','PEG'}` | named solvent/cryoprotectant only | keeping ligands and metals |
| By blanket HETFLAG | `r.id[0] != ' '` | ALL hetero incl. Zn/Mg/heme/FAD/MSE | almost never - breaks binding sites |
| By residue name (naive) | `r.resname == 'HOH'` | misses `'W'`-flagged waters, keeps some | avoid - HETFLAG is authoritative |

## Transforming Coordinates

```python
from Bio.PDB import PDBParser, PDBIO
import numpy as np

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'protein.pdb')

# Entity.transform applies coords @ rot + tran (row convention) to every atom in place.
identity = np.identity(3)
translation = np.array([10.0, 0.0, 0.0])
structure.transform(identity, translation)

io = PDBIO()
io.set_structure(structure)
io.save('translated.pdb')
```

## Rotation Around an Axis

```python
from Bio.PDB import PDBParser
from Bio.PDB.vectors import rotaxis, Vector
import numpy as np

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'protein.pdb')

# rotaxis returns a row-convention matrix intended for Entity/atom.transform.
rot = rotaxis(np.radians(90), Vector(0, 0, 1))

# Rotate about the center of mass: pick tran so the center is the fixed point of coords @ rot + tran.
center = np.array([a.coord for a in structure.get_atoms()]).mean(axis=0)
tran = center - center @ rot
structure.transform(rot, tran)
```

## Applying an External / Assembly Operator

```python
from Bio.PDB import PDBParser
import numpy as np

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'protein.pdb')

# REMARK 350 / _pdbx_struct_oper_list operators are column-convention: newcoord = R @ coord + t.
R = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
t = np.array([25.0, 0.0, 0.0])

# Entity.transform expects the row convention, so transpose the column-convention R first.
structure.transform(R.T, t)
```

## Center Structure at Origin

```python
from Bio.PDB import PDBParser
import numpy as np

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'protein.pdb')

center = np.array([a.coord for a in structure.get_atoms()]).mean(axis=0)
structure.transform(np.identity(3), -center)
```

## Removing Atoms, Residues, and Chains

```python
from Bio.PDB import PDBParser, PDBIO

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'protein.pdb')
model = structure[0]

# Detach hydrogens; collect ids first so the child dict is not mutated mid-iteration.
for residue in model.get_residues():
    for atom_id in [a.id for a in residue if a.element == 'H']:
        residue.detach_child(atom_id)

# Detach whole chains by id.
if model.has_id('B'):
    model.detach_child('B')

io = PDBIO()
io.set_structure(structure)
io.save('cleaned.pdb')
```

## Stripping Solvent by HETFLAG

**Goal:** Remove crystallographic water without deleting functional heteroatoms.

**Approach:** Filter on the residue-id HETFLAG (`r.id[0]`), which is `'W'` for water and `'H_<name>'` for other hetero groups - not on the residue name, which silently keeps `'W'`-flagged waters and cannot distinguish a catalytic metal from a buffer ion.

```python
from Bio.PDB import PDBParser, PDBIO

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'protein.pdb')

# 'W' HETFLAG isolates water; a blanket r.id[0] != ' ' would also delete Zn/Mg/heme/FAD and MSE.
for chain in structure[0]:
    for res_id in [r.id for r in chain if r.id[0] == 'W']:
        chain.detach_child(res_id)

io = PDBIO()
io.set_structure(structure)
io.save('no_water.pdb')
```

## Extracting a Selection with PDBIO Select

```python
from Bio.PDB import PDBParser, PDBIO, Select

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'protein.pdb')

# Select writes a filtered copy without mutating the parsed tree.
class CoreChain(Select):
    def accept_chain(self, chain):
        return chain.id == 'A'
    def accept_residue(self, residue):
        return residue.id[0] == ' ' and 50 <= residue.id[1] <= 100

io = PDBIO()
io.set_structure(structure)
io.save('coreA_50_100.pdb', CoreChain())
```

## Overloading the B-factor Column (Destructive)

**Goal:** Paint a per-residue scalar (conservation, pLDDT) into the B-factor column for viewer coloring.

**Approach:** Overwriting `atom.bfactor` DESTROYS the real temperature factors (and for AlphaFold models overwrites the pLDDT already stored there), so snapshot the originals before writing, set the score on EVERY atom of the residue, and let the viewer autoscale rather than hand-scaling.

```python
from Bio.PDB import PDBParser, PDBIO

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'protein.pdb')

# Snapshot originals: this column is a real temperature factor (or AlphaFold pLDDT) until overwritten.
original_bfactors = {atom.get_full_id(): atom.bfactor for atom in structure.get_atoms()}

conservation = {100: 9.0, 101: 5.0, 102: 3.0}
for residue in structure.get_residues():
    score = conservation.get(residue.id[1])
    if score is None:
        continue
    for atom in residue:
        atom.bfactor = score  # set on all atoms so per-atom coloring is not patchy

io = PDBIO()
io.set_structure(structure)
io.save('colored.pdb')  # do not feed this file back to refinement/validation
```

## Modifying Occupancy

```python
from Bio.PDB import PDBParser, PDBIO

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'protein.pdb')

# Occupancy must stay consistent with altlocs: complementary altlocs should sum to <= 1.
for atom in structure[0]['A'].get_atoms():
    atom.occupancy = 1.0

io = PDBIO()
io.set_structure(structure)
io.save('occupancy_set.pdb')
```

## Renumbering Residues

A sequential renumber like the one below is safe ONLY for internal bookkeeping. To renumber a structure so it matches the UniProt CANONICAL numbering (for figures or mutation mapping), a sequential or fixed-offset renumber SILENTLY MISALIGNS wherever the construct has an expression tag, an unresolved N-terminus, an engineered mutation, or a missing-density loop - which is almost always. Map residue-by-residue through SIFTS / the author `auth_seq_id` scheme instead (see structure-navigation for the observed-vs-SEQRES-vs-UniProt distinction and database-access/uniprot-access for the SIFTS mapping); never assume position N in the file is UniProt residue N.

```python
from Bio.PDB import PDBParser, PDBIO

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'protein.pdb')
chain = structure[0]['A']

# Preserve the (hetflag, ..., icode) tuple; only the resseq middle field changes.
# Assign into a temporary range first to avoid colliding with existing ids mid-loop.
for offset, residue in enumerate(list(chain)):
    hetflag, _, icode = residue.id
    residue.id = (hetflag, offset + 10000, icode)
for new_seq, residue in enumerate(list(chain), start=1):
    hetflag, _, icode = residue.id
    residue.id = (hetflag, new_seq, icode)

io = PDBIO()
io.set_structure(structure)
io.save('renumbered.pdb')
```

## Building a Structure with StructureBuilder

**Goal:** Construct a valid SMCRA tree from coordinates alone.

**Approach:** `StructureBuilder` wires the Structure > Model > Chain > Residue > Atom parent-child links automatically, which is why the writer emits valid records - hand-assembling `Atom` objects without `add` leaves orphans.

```python
from Bio.PDB import StructureBuilder, PDBIO
import numpy as np

sb = StructureBuilder.StructureBuilder()
sb.init_structure('built')
sb.init_model(0)
sb.init_chain('A')
sb.init_seg(' ')
sb.init_residue('ALA', ' ', 1, ' ')
sb.init_atom('N', np.array([-1.0, 0.0, 0.0]), 20.0, 1.0, ' ', 'N', 1, 'N')
sb.init_atom('CA', np.array([0.0, 0.0, 0.0]), 20.0, 1.0, ' ', 'CA', 2, 'C')
sb.init_atom('C', np.array([1.0, 0.0, 0.0]), 20.0, 1.0, ' ', 'C', 3, 'C')
sb.init_atom('O', np.array([1.5, 1.0, 0.0]), 20.0, 1.0, ' ', 'O', 4, 'O')

io = PDBIO()
io.set_structure(sb.get_structure())
io.save('built_structure.pdb')
```

## Copying a Chain (Preserving SMCRA Links)

```python
from Bio.PDB import PDBParser, PDBIO
import copy

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'protein.pdb')

# deepcopy carries the whole subtree with intact parent-child links; reassign id and detach the old parent.
new_chain = copy.deepcopy(structure[0]['A'])
new_chain.id = 'B'
new_chain.detach_parent()
structure[0].add(new_chain)

io = PDBIO()
io.set_structure(structure)
io.save('duplicated_chain.pdb')
```

## Merging Two Structures Without ID Collisions

```python
from Bio.PDB import PDBParser, PDBIO
import copy

parser = PDBParser(QUIET=True)
struct1 = parser.get_structure('s1', 'structure1.pdb')
struct2 = parser.get_structure('s2', 'structure2.pdb')

# Assign explicit non-colliding ids from a free pool; chr(ord(id)+10) breaks on multi-char/adjacent ids.
used = {c.id for c in struct1[0]}
free = (c for c in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' if c not in used)
for chain in list(struct2[0]):
    moved = copy.deepcopy(chain)
    moved.id = next(free)
    moved.detach_parent()
    struct1[0].add(moved)

io = PDBIO()
io.set_structure(struct1)
io.save('merged.pdb')
```

## Common Errors

| Symptom | Cause | Fix |
|---|---|---|
| Rotated structure looks mirrored or points the wrong way | Column-convention operator (REMARK 350 / `_pdbx_struct_oper_list`) applied with the row-convention `Entity.transform` | Transpose first: `structure.transform(R.T, t)`; or apply `R @ coord + t` explicitly |
| Superimposer rotation gives garbage when reused via `np.dot(rot, coord)` | `Superimposer.rotran` is row-convention (`coords @ rot`); `np.dot(rot, coord)` applies the transpose | Use `Entity.transform(rot, tran)` or `coord @ rot + tran` |
| Original structure changed after a transform | All edits mutate in place; a reference is not a copy | `copy.deepcopy(structure)` before modifying |
| B-factors lost / AlphaFold confidence gone after coloring | Writing a scalar into `atom.bfactor` overwrites the temperature factor (or pLDDT) | Snapshot originals first; never send the overloaded file to refinement |
| Catalytic metal or cofactor missing after "removing hetero" | Stripped by `r.id[0] != ' '` or by residue name, deleting Zn/Mg/heme/MSE | Strip water only (`r.id[0] == 'W'`) or use an explicit deny-list |
| `RuntimeError: dictionary changed size during iteration` | Detaching children while iterating the parent | Collect ids into a list first, then `detach_child` |
| `KeyError` when accessing a renumbered residue | Reduced id to `id[1]`, dropping the `(hetflag, ..., icode)` tuple | Key on the full tuple; only display `id[1]` |
| Writer emits truncated or duplicate records | Renumber/merge produced a colliding `(hetflag, resseq, icode)` or chain id | Renumber via a temporary offset; assign ids from a checked free pool |
| Built structure writes an empty or broken file | `Atom`/`Residue` objects created without `add`, leaving SMCRA links unset | Use `StructureBuilder` or wire `add` at every level |
| Only one alternate conformer written after occupancy edit | Altloc/occupancy edited independently so occupancies no longer sum to <= 1 | Keep complementary altlocs consistent as a pair |
| Chain-merge crashes on multi-character chain ids | `chr(ord(chain.id) + 10)` assumes single adjacent characters | Assign explicit ids from a free-id pool |
| mmCIF metadata or anisotropic B-factors dropped after a Bio.PDB round-trip | Bio.PDB does not round-trip ANISOU or the full mmCIF model | For mmCIF-fidelity edits use gemmi; keep Bio.PDB for PDB-scale work |

## Related Skills

- structure-io - Parse and write structure files; mmCIF vs PDB format ceilings
- structure-navigation - Walk chains/residues/atoms and the SMCRA id tuple; observed-vs-SEQRES-vs-UniProt numbering before renumbering
- database-access/uniprot-access - SIFTS mapping of structure residues to UniProt canonical numbering (do not renumber sequentially)
- geometric-analysis - Superimpose structures and read back the row-convention rotation
- interface-analysis - Analyze interfaces after generating the biological assembly
- structure-preparation - Add hydrogens, protonation states, and missing atoms (this skill only removes/edits)
- sequence-manipulation/seq-objects - Generate sequences from modified structures

## References

- Hamelryck T, Manderick B. 2003. PDB file parser and structure class implemented in Python. *Bioinformatics* 19(17):2308-2310. doi:10.1093/bioinformatics/btg332
- Cock PJA, Antao T, Chang JT, et al. 2009. Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics* 25(11):1422-1423. doi:10.1093/bioinformatics/btp163
- Berman HM, Westbrook J, Feng Z, et al. 2000. The Protein Data Bank. *Nucleic Acids Res* 28(1):235-242. doi:10.1093/nar/28.1.235
- wwPDB / RCSB PDB. Biological assembly operators (REMARK 350; `_pdbx_struct_assembly_gen` and `_pdbx_struct_oper_list`). https://www.rcsb.org/docs/programmatic-access/file-download-services

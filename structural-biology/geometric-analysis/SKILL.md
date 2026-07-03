---
name: bio-structural-biology-geometric-analysis
description: Measures geometric properties of protein structures with Biopython Bio.PDB - interatomic distances, distance matrices, bond and dihedral angles (phi/psi/chi, Ramachandran), superposition and RMSD, center of mass, radius of gyration, and solvent accessible surface area (SASA). Use when deciding that RMSD depends on BOTH the superposition and the atom selection (a global all-atom RMSD is dominated by flexible loops and hinge motion and is NOT a cross-protein similarity metric); choosing the metric that matches the question (RMSD for same-molecule displacement, TM-score for same-fold, lDDT for superposition-free local model quality - the quantity pLDDT predicts); recognizing Superimposer needs an equal-length ordered atom-to-atom correspondence; and reporting SASA only alongside its probe radius (1.4A water, Shrake-Rupley) with a preference for relative SASA. Keywords RMSD, TM-score, lDDT, SASA, Shrake-Rupley, superposition, Kabsch, dihedral, Ramachandran, radius of gyration.
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

# Geometric Analysis

**"Calculate the RMSD between two conformations"** -> superimpose on a defined atom correspondence, then report deviation.
- Python: `Bio.PDB.Superimposer` (SVD/Kabsch); `Bio.PDB.qcprot` for speed in tight loops
**"How buried is this residue?"** -> compute solvent accessible surface area and normalize to a per-residue maximum.
- Python: `Bio.PDB.SASA.ShrakeRupley`, then relative SASA against a max-ASA scale

## Governing Principle: the metric IS the question

RMSD is NOT a property of two structures. It is a property of two structures GIVEN a superposition AND an atom selection - change either and the number changes. Report what was aligned (CA-only? a defined core? all-atom?) or the number is uninterpretable.

Global all-atom RMSD is a mean of SQUARED per-atom deviations after a least-squares rigid-body fit, so it is dominated by the worst-fitting atoms. A structure whose 150-residue core is essentially identical but whose 10-residue loop or a hinge-rotated domain swings out by 15A reports a "bad" whole-molecule RMSD (often 4-8A) that hides a near-perfect core. RMSD is also length-dependent (longer proteins accumulate larger RMSD for the same local quality) and is NOT a cross-protein similarity metric - it is only meaningful when a genuine 1:1 correspondence exists (same protein, two states; a model vs its native).

`Superimposer` requires an equal-length, ORDERED atom-to-atom correspondence. Feeding it mismatched or unequal atom lists is the classic error; it computes the optimal fit (Kabsch via SVD), it does NOT solve which atom maps to which. Structures with different sequences need a structure-based alignment FIRST to establish the correspondence (see alignment/structural-alignment), then superposition.

SASA depends on the PROBE RADIUS (1.4A water is a convention, not a constant of nature), the algorithm (Bio.PDB `ShrakeRupley` is Shrake-Rupley only; `freesasa` offers Lee-Richards and LCPO), and whether hydrogens are present. A SASA number without its probe radius is meaningless, and absolute SASA in A^2 is not portable across tools. Prefer RELATIVE SASA (residue SASA / max-ASA of that residue type) using the Tien et al 2013 max-ASA scale.

Backbone phi/psi and omega/cis-peptides are a VALIDATION signal, not a description: a residue in a sterically disallowed Ramachandran region usually means a modeling error, not exotic biology. This skill computes the angles; interpreting outliers as quality flags belongs to structural-biology/structure-validation.

### Decision: which comparison metric

| Metric | Answers (use for) | Caveat | Who reports it |
|---|---|---|---|
| RMSD on a defined core | same molecule, how far did it move after best-fit | needs a real 1:1 correspondence; outlier-dominated (squared mean); length- and selection-dependent; NOT cross-protein | Bio.PDB `Superimposer.rms` |
| TM-score (>0.5 = same fold) | different proteins - same fold? fold recognition | length-normalized and outlier-resistant, but asymmetric (state the reference chain); needs an alignment first | TM-align / US-align (alignment/structural-alignment) |
| GDT-TS / GDT-HA | CASP-style full-model accuracy vs native | superposition-based but fraction-within-cutoff, not a squared mean | LGA / CASP assessors |
| lDDT (0-100; pLDDT for AlphaFold models) | local model quality, multi-domain, without picking a superposition | superposition-free, immune to domain motion; the quantity pLDDT predicts | OpenStructure lDDT / AlphaFold pLDDT |

For cross-protein or fold-similarity work, do not stretch RMSD - route to alignment/structural-alignment (TM-align / Foldseek / DALI). For Ramachandran/omega as a quality gate, route to structural-biology/structure-validation.

## Distance Between Atoms

```python
from Bio.PDB import PDBParser
import numpy as np

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'protein.pdb')

chain = structure[0]['A']
atom1, atom2 = chain[100]['CA'], chain[200]['CA']

distance = atom1 - atom2                              # Atom subtraction returns the distance directly
print(f'Distance: {distance:.2f} A')
print(np.linalg.norm(atom1.coord - atom2.coord))     # Equivalent via numpy on the .coord arrays
```

## Distance Matrix

```python
import numpy as np
from Bio.PDB import PDBParser

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'protein.pdb')

ca_atoms = [r['CA'] for r in structure.get_residues() if r.has_id('CA') and r.id[0] == ' ']  # id[0]==' ' drops waters/hetero
n = len(ca_atoms)

dist = np.zeros((n, n))
for i in range(n):
    for j in range(i + 1, n):
        dist[i, j] = dist[j, i] = ca_atoms[i] - ca_atoms[j]
print(f'Distance matrix: {dist.shape}')
```

## Bond Angle

```python
import numpy as np
from Bio.PDB import PDBParser, calc_angle

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'protein.pdb')

res = structure[0]['A'][100]
angle = calc_angle(res['N'].get_vector(), res['CA'].get_vector(), res['C'].get_vector())  # calc_angle needs Vector, not .coord
print(f'N-CA-C angle: {np.degrees(angle):.1f} deg')
```

## Backbone Dihedrals (phi / psi)

```python
import numpy as np
from Bio.PDB import PDBParser, calc_dihedral

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'protein.pdb')
chain = structure[0]['A']

prev, curr, nxt = chain[99], chain[100], chain[101]
phi = calc_dihedral(prev['C'].get_vector(), curr['N'].get_vector(), curr['CA'].get_vector(), curr['C'].get_vector())
psi = calc_dihedral(curr['N'].get_vector(), curr['CA'].get_vector(), curr['C'].get_vector(), nxt['N'].get_vector())
print(f'phi={np.degrees(phi):.1f}  psi={np.degrees(psi):.1f}')
```

## Ramachandran Angles for All Residues

```python
import numpy as np
from Bio.PDB import PDBParser, PPBuilder

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'protein.pdb')

ppb = PPBuilder()                                    # PPBuilder builds peptides from connectivity, so chain breaks split them
rama = []
for pp in ppb.build_peptides(structure):
    for res, (phi, psi) in zip(pp, pp.get_phi_psi_list()):
        if phi is not None and psi is not None:      # None at termini and chain breaks by design - skip, do not fabricate
            rama.append((res.resname, np.degrees(phi), np.degrees(psi)))
print(f'{len(rama)} residues with phi/psi')
```

Outliers in disallowed regions are usually refinement errors, not biology - interpret them as a quality gate in structural-biology/structure-validation.

## Chi Angles (Sidechain Dihedrals)

```python
import numpy as np
from Bio.PDB import PDBParser, calc_dihedral

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'protein.pdb')
res = structure[0]['A'][100]

if res.has_id('CB') and res.has_id('CG'):            # Gly/Ala lack CB/CG; chi atom quartets are residue-type specific
    chi1 = calc_dihedral(res['N'].get_vector(), res['CA'].get_vector(), res['CB'].get_vector(), res['CG'].get_vector())
    print(f'Chi1: {np.degrees(chi1):.1f} deg')
```

## Superimposing Structures and RMSD

```python
from Bio.PDB import PDBParser, Superimposer

parser = PDBParser(QUIET=True)
ref = parser.get_structure('ref', 'reference.pdb')
mob = parser.get_structure('mobile', 'mobile.pdb')

ref_ca = [r['CA'] for r in ref.get_residues() if r.has_id('CA') and r.id[0] == ' ']
mob_ca = [r['CA'] for r in mob.get_residues() if r.has_id('CA') and r.id[0] == ' ']
n = min(len(ref_ca), len(mob_ca))                    # Superimposer needs EQUAL-LENGTH ORDERED lists; it does NOT solve correspondence
ref_ca, mob_ca = ref_ca[:n], mob_ca[:n]              # Naive truncation is only valid when residues already correspond 1:1

sup = Superimposer()
sup.set_atoms(ref_ca, mob_ca)                        # Optimal rigid-body fit via SVD/Kabsch
print(f'RMSD (CA): {sup.rms:.2f} A')
rotation, translation = sup.rotran                   # The fitted transform, for reuse on other atoms
sup.apply(mob.get_atoms())                           # Mutates mob in place - copy first if the originals are still needed
```

QCP alternative for speed in tight loops (MD, all-vs-all): `from Bio.PDB.qcprot import QCPSuperimposer` (module `Bio.PDB.qcprot`; historically `Bio.PDB.QCPSuperimposer`), same `set_atoms` / `.rms` / `.rotran` / `apply` interface and identical optimum.

## Per-Residue Deviation After Superposition

Fitting minimizes the squared mean, so a global scalar hides where the structures actually differ. A per-residue deviation plot exposes the outlier domination directly.

```python
import numpy as np
from Bio.PDB import PDBParser, Superimposer

parser = PDBParser(QUIET=True)
ref = parser.get_structure('ref', 'reference.pdb')
mob = parser.get_structure('mobile', 'mobile.pdb')

ref_ca = [r['CA'] for r in ref.get_residues() if r.has_id('CA') and r.id[0] == ' ']
mob_ca = [r['CA'] for r in mob.get_residues() if r.has_id('CA') and r.id[0] == ' ']
n = min(len(ref_ca), len(mob_ca))
ref_ca, mob_ca = ref_ca[:n], mob_ca[:n]

sup = Superimposer()
sup.set_atoms(ref_ca, mob_ca)
sup.apply([a for a in mob_ca])                        # Move only the paired CA set into the fitted frame
deviation = np.array([r - m for r, m in zip(ref_ca, mob_ca)])
print(f'core (<2A) residues: {(deviation < 2.0).sum()} / {n}')   # 2A is a common rigid-core cutoff, not a law
print(f'max deviation: {deviation.max():.2f} A at index {deviation.argmax()}')
```

## Center of Mass and Radius of Gyration

```python
import numpy as np
from Bio.PDB import PDBParser

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'protein.pdb')

atoms = list(structure.get_atoms())
coords = np.array([a.coord for a in atoms])
masses = np.array([{'C': 12.0, 'N': 14.0, 'O': 16.0, 'S': 32.0, 'H': 1.0}.get(a.element, 12.0) for a in atoms])

com = (masses[:, None] * coords).sum(axis=0) / masses.sum()
print(f'Center of mass: {com}')

rg = np.sqrt(np.mean(np.sum((coords - coords.mean(axis=0)) ** 2, axis=1)))  # Unweighted radius of gyration
print(f'Radius of gyration: {rg:.2f} A')
```

## Vector Operations

```python
from Bio.PDB import PDBParser

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'protein.pdb')

v1 = structure[0]['A'][100]['CA'].get_vector()
v2 = structure[0]['A'][101]['CA'].get_vector()

diff = v2 - v1
print(f'length: {diff.norm():.2f}  unit: {diff.normalized()}')
cross = v1 ** v2                                      # ** is cross product on Vector objects
dot = v1 * v2                                         # * is dot product on Vector objects
```

## Solvent Accessible Surface Area (SASA)

```python
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'protein.pdb')

sr = ShrakeRupley(probe_radius=1.40, n_points=100)   # 1.4A = water radius (convention); a SASA number is meaningless without its probe radius
sr.compute(structure, level='R')                     # level R attaches .sasa on each residue; children sum to parents
print(f'total SASA: {sum(r.sasa for r in structure.get_residues() if hasattr(r, "sasa")):.1f} A^2')
```

## Relative SASA and Burial

Absolute SASA is not portable across tools. For burial, normalize to a per-residue maximum (Tien et al 2013 theoretical Gly-X-Gly max-ASA).

```python
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

MAX_ASA = {'ALA': 129.0, 'ARG': 274.0, 'ASN': 195.0, 'ASP': 193.0, 'CYS': 167.0, 'GLU': 223.0, 'GLN': 225.0, 'GLY': 104.0, 'HIS': 224.0, 'ILE': 197.0, 'LEU': 201.0, 'LYS': 236.0, 'MET': 224.0, 'PHE': 240.0, 'PRO': 159.0, 'SER': 155.0, 'THR': 172.0, 'TRP': 285.0, 'TYR': 263.0, 'VAL': 174.0}

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'protein.pdb')
ShrakeRupley().compute(structure, level='R')

buried = 0
for res in structure.get_residues():
    if res.resname in MAX_ASA and hasattr(res, 'sasa'):
        rsa = res.sasa / MAX_ASA[res.resname]
        if rsa < 0.20:                               # RSA < 0.20 is the common buried heuristic (a rule of thumb, not a law)
            buried += 1
print(f'buried residues (RSA < 0.20): {buried}')
```

## Secondary-Structure Assignment (DSSP)

Secondary structure is an INTERPRETATION, not a value stored in the file: DSSP, STRIDE, and P-SEA legitimately disagree by 1-2 residues at helix and strand termini, so name the tool and version and never mix assignments from two tools in one analysis. DSSP places the backbone amide hydrogen itself and scores an electrostatic H-bond energy, so it needs no explicit hydrogens, and it processes only the FIRST model of an ensemble. The binary was renamed `dssp` -> `mkdssp` (v4) and must be installed separately.

```python
from Bio.PDB import PDBParser, DSSP

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'protein.pdb')
model = structure[0]                                  # DSSP runs on ONE model only

dssp = DSSP(model, 'protein.pdb', dssp='mkdssp')      # pass the current binary name explicitly
codes = [dssp[k][2] for k in dssp.keys()]             # 8-state: H G I helix, E B strand, T S P - other
helix = sum(c in 'HGI' for c in codes)
strand = sum(c in 'EB' for c in codes)
print(f'helix {helix}, strand {strand}, other {len(codes) - helix - strand} of {len(codes)}')
```

## Common Errors

| Symptom | Cause | Fix |
|---|---|---|
| `Superimposer` errors on differing list sizes | atom lists unequal length; it needs a 1:1 ordered correspondence | match residues by id, or for sequence-different structures align first (alignment/structural-alignment) |
| DSSP helix/strand counts differ from another tool | DSSP, STRIDE, P-SEA disagree at element termini; no ground truth | name the tool+version; never mix assignments; compare like-for-like |
| RMSD is 4-8A for structures that clearly share a fold | global fit dominated by flexible loops/termini/hinge; mean of SQUARED deviations | fit on a defined rigid core, report core vs mobile separately; or use per-residue deviation / TM-score |
| RMSD differs between runs on the "same" pair | different atom selection (CA vs all-atom) or superposition | state the correspondence and fit selection explicitly and hold it constant |
| Ranking models of different length by RMSD | RMSD is length-dependent and not a cross-protein metric | use TM-score (length-normalized) or lDDT (superposition-free) |
| `calc_angle`/`calc_dihedral` AttributeError | passed numpy arrays (`.coord`) not `Vector` objects | pass `atom.get_vector()` |
| phi/psi is `None` at chain ends | terminal residues lack a preceding C or following N | skip `None`; `get_phi_psi_list` returns `None` at breaks/termini by design |
| SASA disagrees with a published value | different probe radius, radii set, algorithm, or H atoms present | recompute all structures like-for-like in one tool; report probe_radius; prefer relative SASA |
| `.sasa` attribute missing on residues | `compute()` run at the wrong level or read before it | call `sr.compute(entity, level='R')` then read `residue.sasa` |
| Distance matrix polluted by waters/heteroatoms | iterating residues without filtering the hetflag | filter `residue.id[0] == ' '` |
| Chi1 computed for Gly/Ala | those residues have no CB/CG | guard `has_id('CB') and has_id('CG')`; chi quartets are residue-type specific |
| Two crystal forms called "different states" at 2A RMSD | difference within coordinate uncertainty / ensemble spread | compare against B-factors, resolution, and NMR ensemble spread before claiming a state change |
| Calling a predicted model "wrong" where it deviates from a crystal structure | the deviating region may be low-pLDDT, a PAE-uncertain inter-domain float, or the crystal is a different (holo/packing) state | overlay pLDDT/PAE on the deviation before judging (alphafold-predictions); confirm it is not just a state difference |
| Original coordinates changed unexpectedly | `Superimposer.apply` and `atom.transform` mutate in place | copy the structure first if the untransformed coordinates are still needed |

## Related Skills

- structure-io - Parse and write PDB/mmCIF structure files
- structure-navigation - Walk chains, residues, atoms; handle altlocs and disordered residues
- structure-modification - Transform coordinates and edit structures in place
- structural-biology/interface-analysis - Residue contacts, contact maps, and buried-surface interface analysis (NeighborSearch)
- structural-biology/structure-validation - Ramachandran and omega/cis-peptide outliers as a quality gate
- structural-biology/alphafold-predictions - overlay pLDDT/PAE when a compared structure is a predicted model
- structural-biology/modern-structure-prediction - reconcile predicted models via pLDDT/PAE/pTM before RMSD claims
- alignment/structural-alignment - Cross-protein fold comparison and correspondence (TM-align, Foldseek, DALI)

## References

- Cock PJA, et al. (2009) Biopython. *Bioinformatics* 25(11):1422-1423.
- Kabsch W (1976) A solution for the best rotation to relate two sets of vectors. *Acta Crystallogr A* 32:922-923.
- Theobald DL (2005) Rapid calculation of RMSDs using a quaternion-based characteristic polynomial. *Acta Crystallogr A* 61(4):478-480.
- Zhang Y, Skolnick J (2004) Scoring function for automated assessment of protein structure template quality. *Proteins* 57(4):702-710.
- Xu J, Zhang Y (2010) How significant is a protein structure similarity with TM-score = 0.5? *Bioinformatics* 26(7):889-895.
- Mariani V, Biasini M, Barbato A, Schwede T (2013) lDDT: a local superposition-free score for comparing protein structures and models. *Bioinformatics* 29(21):2722-2728.
- Shrake A, Rupley JA (1973) Environment and exposure to solvent of protein atoms. Lysozyme and insulin. *J Mol Biol* 79(2):351-371.
- Tien MZ, Meyer AG, Sydykova DK, Spielman SJ, Wilke CO (2013) Maximum allowed solvent accessibilities of residues in proteins. *PLoS ONE* 8(11):e80635.

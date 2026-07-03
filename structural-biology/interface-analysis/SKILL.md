---
name: bio-structural-biology-interface-analysis
description: Maps protein-protein and protein-ligand interfaces with Bio.PDB, computing contact residues and buried surface area (BSA). Use when choosing a contact cutoff and stating its rationale (heavy-atom 4-5A vs CA-CA 8A vs a SASA-based definition); deciding a contact list is not an interface and computing buried surface area (dSASA/BSA) instead; distinguishing a genuine biological interface from a crystal-packing artifact; identifying ligand-contact or epitope residues; and computing on the biological assembly rather than the asymmetric unit. Keywords interface, buried surface area, BSA, contacts, NeighborSearch, PISA, crystal packing, epitope, binding site, ShrakeRupley.
tool_type: python
primary_tool: Bio.PDB
---

## Version Compatibility

Reference examples tested with: biopython 1.83+, numpy 1.26+, freesasa 2.2+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Interface Analysis

**"Which residues contact the ligand / the partner chain?"** -> Threshold interatomic distances and collect residues within a cutoff.
- Python: `Bio.PDB.NeighborSearch(atoms).search_all(cutoff, level='R')` or `.search(center, cutoff, level='R')`

**"Compute the buried surface area of this interface"** -> Subtract complex SASA from the summed SASA of the isolated partners.
- Python: `Bio.PDB.SASA.ShrakeRupley().compute(entity)`, then `BSA = SASA_A + SASA_B - SASA_complex`

**"Is this interface biological or a crystal-packing artifact?"** -> Score interface size and chemistry on the biological assembly, and treat the assignment as a hypothesis.
- Python: Bio.PDB for BSA / H-bonds; PDBePISA for the assembly call (external service)

## Governing Principle

A "contact" is not a physical fact, it is a thresholded distance, and the residue count changes with the cutoff. Heavy-atom pairs within 4-5A capture direct van der Waals contact; CA-CA within 8A captures topological proximity (contact maps, coevolution features) but says nothing about side-chain interaction; 3.5-4.0A heavy-atom is the H-bond / salt-bridge regime. Changing 4A to 5A can shift the contact count substantially, so a contact result is meaningless without its atom set and cutoff stated explicitly (Chakrabarti & Janin 2002 *Proteins* 47:334-343). The trap is reporting "N interface residues" or "N contacts" as if the number were intrinsic.

A contact list is not an interface. The physical interface measure is BURIED SURFACE AREA (BSA, also dSASA): BSA = SASA(part A alone) + SASA(part B alone) - SASA(complex), conventionally halved to report the area buried per partner. All three SASA terms must be computed with identical parameters (same probe radius, radii set, algorithm) or the subtraction is garbage (see geometric-analysis for SASA fundamentals). SASA itself depends on the probe radius (1.4A water default) and the algorithm, so an absolute BSA is only comparable to another BSA computed the same way. Because heavy-atom contacts and BSA are both defined on non-hydrogen atoms, hydrogens are NOT required for either - add them (structure-preparation) only for H-bond/salt-bridge angle geometry, and if H are present keep them consistent across all three SASA terms.

The deepest trap: an interface seen in the deposited ASYMMETRIC UNIT may be a CRYSTAL-PACKING ARTIFACT, not biology. The asymmetric unit is a crystallographic bookkeeping object; the functional molecule is the BIOLOGICAL ASSEMBLY, which may be a subset of the ASU or built from several ASUs by symmetry. Compute interfaces on the biological assembly, not blindly on the ASU (see structure-io for downloading the assembly). PDBePISA (Krissinel & Henrick 2007 *J Mol Biol* 372:774-797) predicts the biological assembly and scores interface stability, but it recovers the correct assembly only ~80-90% of the time and has known false positives, so "biological interface" is a HYPOTHESIS. Larger BSA, more H-bonds and salt bridges, shape complementarity, and evolutionary conservation of interface residues each raise confidence, but each is probabilistic, not proof (Levy 2010 *J Mol Biol* 403:660-670). Corroborate anything load-bearing with solution data (SEC-MALS, SAXS, native MS).

## Decision: contact / interface definition

| Definition | What it captures | Best when | Fails / misleads when |
|---|---|---|---|
| Heavy-atom (non-H) <= 4-5A | Direct physical / vdW contact | Interface residue lists, ligand-contact residues, epitopes | Cutoff unstated; H atoms present shift the count |
| CA-CA <= 8A | Topological proximity of backbones | Contact maps, coevolution / ML features, fold fingerprint | Read as "side chains interact" - it does not imply that |
| Heavy-atom 3.5-4.0A + angle | H-bonds / salt bridges | Chemistry of the interface | Definitions are loose and tool-dependent (state exact criteria) |
| BSA / dSASA (SASA-based) | Physical extent of the interface (area) | Quantifying interface size, biological-vs-crystal | Terms computed with mismatched SASA parameters |

The one-line rule: heavy-atom 4-5A answers "who touches"; CA-CA 8A answers "who is near"; BSA answers "how big is the interface". State the atom set and cutoff every time.

## Decision: biological interface vs crystal contact

| Signal | Biological interface tends to | Crystal contact tends to | Caveat |
|---|---|---|---|
| Buried surface area (per side) | Larger, often > ~800-1000 A^2 | Small, often < ~400 A^2 | Wide overlap; not a hard cutoff |
| H-bonds / salt bridges | More, specific | Few, incidental | Definition-dependent counts |
| Shape complementarity | High | Lower | Not diagnostic alone |
| Interface residue conservation | Conserved across homologs | Not conserved | Needs an alignment / ortholog set |
| PDBePISA assignment | Called stable (CSS toward 1.0) | Called unstable | ~80-90% accurate; known false positives |
| Recurs across crystal forms | Yes | No (packing-specific) | Requires multiple depositions |

Every row is probabilistic. Interface size (BSA) is the single most-used signal, but small biological interfaces (transient/weak complexes) and large crystal contacts both exist, so no one number settles it.

## Contact residues between two chains

**Goal:** List the residues of chain A and chain B that form the interface, under an explicit cutoff.

**Approach:** Build one KD-tree over the interface atoms, query all close pairs at residue level, and keep pairs whose two residues belong to different chains. Heavy-atom cutoff 4.5A (midpoint of the 4-5A vdW-contact regime; excludes H so it is robust to whether H atoms were modeled).

```python
from Bio.PDB import PDBParser, NeighborSearch, Selection

parser = PDBParser(QUIET=True)
structure = parser.get_structure('complex', 'complex.pdb')
model = structure[0]

cutoff = 4.5  # heavy-atom contact; 4-5A captures direct vdW contact, state it always
atoms = [a for a in model.get_atoms() if a.element != 'H']
ns = NeighborSearch(atoms)

interface_a, interface_b = set(), set()
for res1, res2 in ns.search_all(cutoff, level='R'):
    c1, c2 = res1.get_parent().id, res2.get_parent().id
    if c1 == 'A' and c2 == 'B':
        interface_a.add(res1); interface_b.add(res2)
    elif c1 == 'B' and c2 == 'A':
        interface_b.add(res1); interface_a.add(res2)

print(f'Chain A interface residues ({cutoff}A): {len(interface_a)}')
print(f'Chain B interface residues ({cutoff}A): {len(interface_b)}')
```

## Ligand-contact (binding-site / epitope) residues

**Goal:** Identify the protein residues lining a bound ligand or the residues an antibody contacts (structural epitope).

**Approach:** Select the ligand atoms (a HETATM group, hetflag starts with 'H_'), search protein atoms within the cutoff of each, collect unique parent residues. The same pattern with two protein chains yields a structural epitope.

```python
from Bio.PDB import PDBParser, NeighborSearch

parser = PDBParser(QUIET=True)
structure = parser.get_structure('complex', 'complex.pdb')
model = structure[0]

ligand_resname = 'ATP'  # target HETATM group
cutoff = 4.5

ligand_atoms = [a for r in model.get_residues() if r.resname == ligand_resname
                for a in r if a.element != 'H']
protein_atoms = [a for a in model.get_atoms()
                 if a.element != 'H' and a.get_parent().id[0] == ' ']
ns = NeighborSearch(protein_atoms)

pocket = set()
for a in ligand_atoms:
    for res in ns.search(a.coord, cutoff, level='R'):
        pocket.add((res.get_parent().id, res.id[1], res.resname))

for chain, num, name in sorted(pocket):
    print(f'{chain} {name}{num}')
```

## Buried surface area (BSA / dSASA)

**Goal:** Quantify the physical size of a two-chain interface as area buried on complex formation.

**Approach:** Compute SASA on the intact complex, then on each chain in isolation (same ShrakeRupley settings), and take BSA = SASA_A + SASA_B - SASA_complex. Halve for per-partner area. Probe radius 1.4A models a water molecule; keep it identical across all three computations or the subtraction is meaningless.

```python
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley

parser = PDBParser(QUIET=True)
sr = ShrakeRupley(probe_radius=1.4)  # 1.4A ~ water; MUST match across all three terms

def chain_sasa(path, keep_chains):
    structure = parser.get_structure('s', path)
    model = structure[0]
    for chain in list(model):
        if chain.id not in keep_chains:
            model.detach_child(chain.id)
    sr.compute(model, level='C')
    return sum(chain.sasa for chain in model)

sasa_complex = chain_sasa('complex.pdb', {'A', 'B'})
sasa_a = chain_sasa('complex.pdb', {'A'})
sasa_b = chain_sasa('complex.pdb', {'B'})

bsa_total = sasa_a + sasa_b - sasa_complex
print(f'Total buried surface area: {bsa_total:.0f} A^2')
print(f'Per partner: {bsa_total / 2:.0f} A^2')  # convention: split half to each side
```

For Lee-Richards SASA or full control of the radii set and probe, use `freesasa` instead of ShrakeRupley (Mitternacht 2016 *F1000Research* 5:189); Bio.PDB provides only Shrake-Rupley. Compute all three terms in the same tool.

## H-bonds and salt bridges (geometric heuristics)

**Goal:** Estimate the specific polar interactions across an interface.

**Approach:** Salt bridge = an acidic side-chain oxygen (Asp/Glu OD/OE) within ~4A of a basic side-chain nitrogen (Arg/Lys/His NZ/NH/NE/ND). These definitions are loose and tool-dependent; state the exact distance (and any angle) used. Without modeled hydrogens, a true H-bond angle cannot be checked, so the distance-only result is an upper bound.

```python
from Bio.PDB import PDBParser, NeighborSearch

parser = PDBParser(QUIET=True)
model = parser.get_structure('c', 'complex.pdb')[0]

acidic = {('ASP', 'OD1'), ('ASP', 'OD2'), ('GLU', 'OE1'), ('GLU', 'OE2')}
basic = {('ARG', 'NH1'), ('ARG', 'NH2'), ('ARG', 'NE'),
         ('LYS', 'NZ'), ('HIS', 'ND1'), ('HIS', 'NE2')}
salt_cutoff = 4.0  # common salt-bridge distance; literature ranges 3.2-5.0A, report the choice

ns = NeighborSearch(list(model.get_atoms()))
bridges = []
for a1, a2 in ns.search_all(salt_cutoff, level='A'):
    k1 = (a1.get_parent().resname, a1.name)
    k2 = (a2.get_parent().resname, a2.name)
    cross = a1.get_parent().get_parent().id != a2.get_parent().get_parent().id
    if cross and ((k1 in acidic and k2 in basic) or (k1 in basic and k2 in acidic)):
        bridges.append((a1.get_parent(), a2.get_parent()))

print(f'Candidate interchain salt bridges (<= {salt_cutoff}A): {len(bridges)}')
```

## PDBePISA for the biological assembly

PDBePISA computes interfaces and predicts the biological assembly from the crystal, reporting interface area, an interface solvation free energy of assembly, the number of H-bonds and salt bridges, and a Complexation Significance Score (CSS, 0-1) ranking each interface by how much it drives assembly. It is a web service (https://www.ebi.ac.uk/pdbe/pisa/) with per-entry results; there is no Bio.PDB binding. Use it to get the assembly call and interface energetics, then treat the assignment as a hypothesis to corroborate (see the biological-vs-crystal table). Do not report the PISA assembly as ground truth.

## Common Errors

| Symptom | Cause | Fix |
|---|---|---|
| Contact count changes between runs / papers | Cutoff or atom set not stated or not matched | Fix and report the cutoff and whether H atoms are included |
| "Interface" that vanishes in solution | Computed on the asymmetric unit, not the biological assembly | Download and compute on the biological assembly (structure-io) |
| BSA comes out near zero or negative | SASA terms computed with different parameters or on different files | Use identical ShrakeRupley settings for complex and each isolated part |
| Huge BSA but no biology | Large crystal contact misread as biological | Cross-check H-bonds, conservation, PISA CSS, recurrence across crystal forms |
| Ligand-contact residues missing | Ligand skipped because it is a HETATM, filtered out with waters | Select the ligand by resname/hetflag before filtering standard residues |
| Doubled / impossible contacts at one residue | Alternate conformations (altloc) both counted | Select one altloc before contact search (structure-modification) |
| Interface residues span a chain gap oddly | Missing/disordered residues modeled as absent | Reconcile against SEQRES; missing loops are disorder, not a real gap |
| H-bond angles cannot be computed | No hydrogens modeled in the file | Report distance-only heuristics as an upper bound, or add H first |
| Salt-bridge count disagrees with another tool | Distance/angle definition differs between tools | State exact criteria; definitions are not standardized |
| CA-CA 8A "interface" implies side-chain contact | 8A is topological proximity, not physical contact | Use heavy-atom 4-5A for physical contact claims |
| SASA / BSA not comparable to a literature value | Different probe radius, radii set, or algorithm | Recompute both like-for-like in one tool |
| PISA assembly taken as fact | PISA is ~80-90% accurate with known false positives | Treat as a hypothesis; corroborate with solution data |

## Related Skills

- geometric-analysis - SASA fundamentals, NeighborSearch, distances that this skill builds on
- structure-io - download the biological assembly (not just the asymmetric unit) before interface analysis
- structure-modification - resolve altlocs and strip waters/additives before contact detection
- structure-navigation - select chains, residues, and HETATM ligands by identity
- structure-validation - check the region of interest is well-fit before trusting an interface
- structure-preparation - add hydrogens before checking H-bond/salt-bridge geometry at an interface
- binding-site-detection - de-novo cavity/pocket detection on apo structures (complement to mapping a bound ligand)
- immunoinformatics/epitope-prediction - structural epitope mapping from antibody-antigen complexes
- chemoinformatics/virtual-screening - binding-site definition for docking
- alignment/structural-alignment - superpose complexes before comparing interfaces

## References

- Krissinel E, Henrick K (2007) Inference of macromolecular assemblies from crystalline state. *J Mol Biol* 372(3):774-797. (PISA / PDBePISA; biological-assembly prediction and its failure modes)
- Chakrabarti P, Janin J (2002) Dissecting protein-protein recognition sites. *Proteins* 47(3):334-343. (interface core/rim dissection; contact-definition dependence)
- Levy ED (2010) A simple definition of structural regions in proteins and its use in analyzing interface evolution. *J Mol Biol* 403(4):660-670. (core-rim-support model; 25% RSA burial threshold)
- Shrake A, Rupley JA (1973) Environment and exposure to solvent of protein atoms. Lysozyme and insulin. *J Mol Biol* 79(2):351-371. (Shrake-Rupley SASA underlying BSA)
- Tien MZ, Meyer AG, Sydykova DK, Spielman SJ, Wilke CO (2013) Maximum allowed solvent accessibilities of residues in proteins. *PLoS ONE* 8(11):e80635. (max-ASA scale for relative burial of interface residues)
- Mitternacht S (2016) FreeSASA: an open source C library for solvent accessible surface area calculations. *F1000Research* 5:189. (Lee-Richards alternative to ShrakeRupley)
- Cock PJA, et al. (2009) Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics* 25(11):1422-1423. (Bio.PDB toolkit)

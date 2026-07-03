---
name: bio-structural-biology-structure-navigation
description: Navigate the Bio.PDB SMCRA hierarchy (Structure-Model-Chain-Residue-Atom) safely, surfacing the heterogeneity it hides by default. Use when deciding how to handle altloc/DisorderedAtom conformers before a distance or RMSD, indexing residues insertion-code-safe with the full (hetflag, resseq, icode) tuple, choosing the ATOM/observed vs SEQRES/canonical vs UniProt sequence, selecting the right Model for an NMR ensemble, filtering waters/hetero/metals correctly, and reconciling auth vs label numbering. Keywords SMCRA, altloc, DisorderedAtom, insertion code, SEQRES, PPBuilder, auth_seq_id.
tool_type: python
primary_tool: Bio.PDB
goal_approach_exempt: true
---

## Version Compatibility

Reference examples tested with: biopython 1.83+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Structure Navigation

**"Walk the chains, residues, and atoms; pull out the sequence"** -> Traverse the Structure-Model-Chain-Residue-Atom (SMCRA) tree, but treat every level as a lossy projection that hides heterogeneity unless asked otherwise.
- Python: `structure[0]['A'][(' ', 100, ' ')]['CA'].coord` for full-tuple direct access

## Governing Principle: SMCRA is a convenient tree that HIDES heterogeneity by default

The SMCRA model (Structure > Model > Chain > Residue > Atom) is a readable in-memory tree, but its defaults quietly collapse the very heterogeneity that changes the answer. Five traps recur, and none of them raises an error - the code runs and returns a plausible-looking number computed on the wrong thing.

1. A `DisorderedAtom` silently forwards every uncaught call to ONE child - the HIGHEST-OCCUPANCY altloc, not literally altloc 'A'. So `get_atoms()`, `.coord`, distances, clashes, and RMSDs all use a single conformer, invisibly, even when the active site is 60/40 disordered. This is the #1 correctness trap. Enumerate with `is_disordered()` and `disordered_get_list()`; never let both altlocs of one atom enter the same geometric calculation.
2. The residue id is a 3-tuple `(hetflag, resseq, icode)`. Naive `residue.id[1]` drops both the hetero flag and the INSERTION CODE, so antibody residues 100, 100A, 100B collapse onto one key; `chain[100]` works only until an insertion code, a hetero residue, or an altloc residue exists at that number, then raises a KeyError that looks like the residue is missing. Key on the full tuple; reduce to `id[1]` only for display.
3. The sequence PPBuilder extracts is the OBSERVED (ATOM-record) sequence with missing-density gaps silently concatenated away - it is NOT the SEQRES/construct/UniProt canonical sequence. A disordered 12-residue loop becomes 12 vanished characters with no marker, so mapping conservation or alignment columns by string position is off-by-many after the first gap. Map through residue NUMBERS (auth_seq_id) or SIFTS, never by string index.
4. NMR and multi-state files have multiple `Model` objects. Iterating chains without first selecting a model conflates conformers; `structure[0]` silently picks one NMR member and over-claims precision. Ask "how many models and why" first.
5. mmCIF carries two numbering schemes: auth (matches the paper, has insertion codes, can be negative/gapped) and label (gapless 1..N, no icodes). `MMCIFParser` defaults to auth; flip `auth_residues=False` and residue 100 becomes a different residue. Pick one scheme and stay in it.

## Decision: which sequence source

These three "sequences" are routinely conflated; each answers a different question and the wrong one silently misindexes everything downstream.

| Source | What it is | Get it via | Use when | Fails when |
|--------|-----------|------------|----------|------------|
| ATOM / observed | Only residues with modeled coordinates; gaps concatenated away | `PPBuilder().build_peptides()` then `pp.get_sequence()` | Per-atom geometry, contacts, extracting exactly what was resolved | Aligning to UniProt/MSA by string position (gaps shift the frame) |
| SEQRES / declared | Full sequence the depositor says is in the crystal, including unresolved residues | `SeqIO.parse(file, 'pdb-seqres')` or `'cif-seqres'` | Knowing the true construct length, locating missing loops | Assuming every SEQRES residue has coordinates (it does not) |
| UniProt / canonical | The reference biological sequence (no tags, no engineered mutations) | `database-access/uniprot-access` + SIFTS residue mapping | Mapping conservation/domains/mutations onto structure positions | Assuming construct == canonical (tags, point mutations, chimeras differ) |

## Decision: residue selection idiom

Filter on the hetflag (`id[0]`), not the residue name, and know exactly what each idiom keeps and drops.

| Goal | Idiom | Keeps / drops correctly? |
|------|-------|--------------------------|
| Standard amino acids only | `r.id[0] == ' '` | Correct; drops water, ligands, and modified residues |
| Water | `r.id[0] == 'W'` | Correct; water hetflag is `'W'`, NOT `'H_'` - a `startswith('H_')` water strip MISSES water |
| Ligands and hetero groups | `r.id[0].startswith('H_')` | Also catches modified residues (MSE, SEP, PTR) mid-chain - not just free ligands |
| Strip hetero blindly | `r.id[0] != ' '` | DANGEROUS - deletes catalytic metals, cofactors, AND selenomethionine (MSE) out of the chain |
| A specific residue type | `r.resname == 'ARG'` | Fine for type queries; never use resname to classify water vs ligand vs standard |

## Required Imports

```python
from Bio.PDB import PDBParser, MMCIFParser, PPBuilder, CaPPBuilder, Selection
from Bio.Data.PDBData import protein_letters_3to1, protein_letters_3to1_extended
```

## Accessing Hierarchy Levels

```python
parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'protein.pdb')

model = structure[0]                       # first model (see NMR caveat below)
chain = model['A']
residue = chain[(' ', 100, ' ')]           # full id tuple - insertion-code and hetero safe
residue_bare = chain[100]                   # convenience path; breaks on icode/hetero/altloc at 100
atom = residue['CA']
```

## Iterating Over Structure

```python
for model in structure:
    for chain in model:
        for residue in chain:
            hetflag, resseq, icode = residue.id   # keep the whole tuple, not resseq alone
            for atom in residue:
                print(f'{chain.id}:{resseq}{icode}:{atom.name}')

for chain in structure.get_chains():
    print(f'Chain: {chain.id}')
```

## Enumerating Disordered Atoms (the highest-value pattern)

```python
# A DisorderedAtom forwards uncaught calls to its highest-occupancy child by default,
# so plain iteration measures ONE conformer. Enumerate every altloc explicitly.
for atom in residue:
    if atom.is_disordered():
        for alt in atom.disordered_get_list():   # each alt is a real Atom with its own coord/occupancy
            print(f'{atom.name} altloc {alt.altloc}: occ={alt.occupancy} coord={alt.coord}')
    else:
        print(f'{atom.name}: coord={atom.coord}')

# disordered_select(altloc) mutates the active child GLOBALLY - reset it when done
if atom.is_disordered():
    atom.disordered_select(atom.disordered_get_id_list()[0])
```

## Disordered Residues (microheterogeneity / point mutation at one site)

```python
# Residue.is_disordered() returns 1 for a NORMAL residue that merely holds altloc atoms and 2 for a
# true DisorderedResidue (two resnames at one position); disordered_get_id_list/disordered_select
# exist ONLY on the latter, so gate on == 2 (or isinstance DisorderedResidue) or this AttributeErrors.
if residue.is_disordered() == 2:
    names = residue.disordered_get_id_list()      # alternative resnames at this position
    residue.disordered_select('ALA')              # pick one by resname before any geometry
```

## Extracting the Observed Sequence (know it has silent gaps)

```python
# PPBuilder builds peptides from OBSERVED atoms via a C-N distance criterion and breaks at gaps;
# the returned Seq has missing-density loops concatenated away with no gap marker.
ppb = PPBuilder()
for pp in ppb.build_peptides(structure):
    print(f'observed segment len={len(pp.get_sequence())}: {pp.get_sequence()}')

# CaPPBuilder connects residues whose CA atoms are within ~4.3 Angstroms, so it bridges
# small backbone breaks - useful for CA-only/broken chains, but it can mis-join true gaps.
ca_ppb = CaPPBuilder()
segments = ca_ppb.build_peptides(structure)
```

## Reading the Declared (SEQRES) Sequence and Locating Gaps

```python
from Bio import SeqIO

# SEQRES = the full declared sequence, including residues with no coordinates.
for record in SeqIO.parse('protein.pdb', 'pdb-seqres'):
    print(f'{record.id} declared length {len(record.seq)}')

# header['missing_residues'] (populated when get_header=True) lists unmodeled residues -
# the difference between SEQRES and observed. Reconcile before mapping to UniProt.
parser = PDBParser(QUIET=True, get_header=True)
structure = parser.get_structure('protein', 'protein.pdb')
missing = structure.header.get('missing_residues', [])
```

## Converting Residue Names (modified residues map to X or drop)

```python
# protein_letters_3to1 is the strict 20-aa map (unknown resname -> KeyError, so use .get).
# protein_letters_3to1_extended additionally maps modified residues (MSE->M, SEP->S, PTR->Y).
seq = ''
for residue in chain:
    if residue.id[0] == ' ' or residue.resname in protein_letters_3to1_extended:
        seq += protein_letters_3to1_extended.get(residue.resname, 'X')
```

## Selecting Entities and Full Identifiers

```python
residues = Selection.unfold_entities(structure, 'R')   # S/M/C/R/A level codes
atoms = Selection.unfold_entities(chain, 'A')

atom = structure[0]['A'][(' ', 100, ' ')]['CA']
print(atom.get_full_id())   # ('protein', 0, 'A', (' ', 100, ' '), ('CA', ' '))
```

## Working with NMR Ensembles (do not conflate models)

```python
n_models = len(structure)                  # NMR = conformer ensemble, X-ray/cryo-EM usually 1
if n_models > 1:
    # Compute per-model and report the distribution; never average coordinates across models.
    for model in structure:
        ca = [r['CA'].coord for r in model.get_residues() if r.has_id('CA')]
        print(f'model {model.id}: {len(ca)} CA atoms')
```

## Reading mmCIF with an Explicit Numbering Scheme

```python
# MMCIFParser defaults to auth numbering (matches the paper, has insertion codes).
# label numbering is gapless 1..N with no icodes - a different residue at the same number.
cif_parser = MMCIFParser(QUIET=True, auth_residues=True)   # keep auth for literature/UniProt cross-ref
structure = cif_parser.get_structure('protein', 'protein.cif')
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Distance/RMSD subtly off on a disordered site | `get_atoms()` returned only the highest-occupancy altloc of a DisorderedAtom | Enumerate `disordered_get_list()`; pick one altloc consistently before geometry |
| Impossibly close contacts / inflated atom count | Both altlocs of one atom entered the same calculation | Select a single altloc per site; never mix conformers |
| KeyError on `chain[100]` for a residue that is clearly present | Residue has an insertion code, hetero flag, or altloc, so the bare-int path misses it | Index with the full tuple `chain[(' ', 100, ' ')]` or iterate and match `id[1]`/`id[2]` |
| Antibody CDR residues 100/100A/100B collapse to one | Keyed on `residue.id[1]` (resseq) and dropped `id[2]` (icode) | Key on the full `(hetflag, resseq, icode)` tuple |
| Structure sequence one residue shorter than expected after each loop | PPBuilder returns the observed sequence with missing-density gaps concatenated away | Use SEQRES (`pdb-seqres`) for length; map to UniProt by residue number or SIFTS, not string index |
| Selenomethionine protein reads full of X or has holes in the chain | MSE is a hetero (`H_MSE`) residue; strict `protein_letters_3to1` returns X or a blanket hetero strip deleted it | Use `protein_letters_3to1_extended` (MSE->M); filter hetero by explicit allow/deny list |
| Water strip leaves waters behind | Filtered with `startswith('H_')`; water hetflag is `'W'` | Strip water with `r.id[0] == 'W'` |
| Stripping hetero removed a catalytic metal or cofactor | Blanket `r.id[0] != ' '` deletes metals, heme, FAD, ions, and mid-chain MSE | Remove only water/buffer by explicit list; keep ligands and metals |
| NMR metric multiplied ~20x or averaged to nonsense | Iterated all models, or averaged coordinates across the ensemble | Select one representative model, or compute per-model and report the spread |
| mmCIF residue numbers do not match the paper | Read label numbering instead of auth (or mixed a label index into an auth structure) | Set `auth_residues=True` (default) and stay in one scheme |
| Missing loop treated as a covalent chain break | Gap in coordinates is unmodeled disorder, not a real break | Reconcile against `header['missing_residues']`/SEQRES; flag gaps as disorder |
| Silent atom drops / merged chains on a messy file | `PDBParser(QUIET=True)` suppressed the discontinuous-chain warnings | For unfamiliar files parse without QUIET (or capture warnings) first |

## Related Skills

- structure-io - Parse and write PDB/mmCIF; auth vs label numbering at the I/O layer
- geometric-analysis - Distances, angles, RMSD, SASA once heterogeneity is resolved
- structure-modification - Strip waters/hetero, edit coordinates and B-factors safely
- interface-analysis - Requires the biological assembly, not the deposited asymmetric unit
- sequence-manipulation/seq-objects - Work with the extracted Seq objects
- alignment/msa-parsing - Map SEQRES/ATOM sequences onto alignment columns
- database-access/uniprot-access - Fetch the canonical sequence for SIFTS-based mapping

## References

- Hamelryck T, Manderick B (2003). PDB file parser and structure class implemented in Python. *Bioinformatics* 19(17):2308-2310. DOI 10.1093/bioinformatics/btg332.
- Cock PJA, Antao T, Chang JT, et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics* 25(11):1422-1423. DOI 10.1093/bioinformatics/btp163.
- Berman HM, Westbrook J, Feng Z, et al. (2000). The Protein Data Bank. *Nucleic Acids Research* 28(1):235-242. DOI 10.1093/nar/28.1.235.
- Velankar S, Dana JM, Jacobsen J, et al. (2013). SIFTS: Structure Integration with Function, Taxonomy and Sequences resource. *Nucleic Acids Research* 41(D1):D483-D489. DOI 10.1093/nar/gks1258.

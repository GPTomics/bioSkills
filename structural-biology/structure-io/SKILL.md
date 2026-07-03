---
name: bio-structural-biology-structure-io
description: Reads, writes, downloads, and converts macromolecular structures with Biopython Bio.PDB. Use when choosing a format (mmCIF/PDBx vs legacy PDB vs BinaryCIF) for a structure that may exceed PDB's ~62-chain / 99,999-atom limits; when residue numbers do not match the paper because of auth_* vs label_* numbering (MMCIFParser defaults auth_residues=True); when metadata (resolution, method, R-free) is missing because Bio.PDB drops it and MMCIF2Dict is needed; when the deposited coordinates are the asymmetric unit and the biological assembly must be downloaded separately; when downloading from RCSB (files.rcsb.org, PDBList); and when a legacy MMTF path is dead (RCSB retired MMTF July 2024, use BinaryCIF).
tool_type: python
primary_tool: Bio.PDB
goal_approach_exempt: true
---

## Version Compatibility

Reference examples tested with: biopython 1.85+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Structure I/O

**"Read a structure file"** -> Parse a deposited coordinate file into an in-memory SMCRA tree, or fetch it from a wwPDB mirror.
- Python: `Bio.PDB.MMCIFParser().get_structure('id', 'file.cif')`, `Bio.PDB.PDBParser()`, `Bio.PDB.PDBList()`

## Governing Principle

mmCIF (PDBx) is the canonical modern format; the legacy fixed-column PDB format is a frozen, lossy container, and a "parse everything as PDB" reflex silently truncates or fails on anything large. The PDB format was frozen in 2012 and cannot physically exceed 99,999 atoms (5-digit serial), ~62 chains (single alphanumeric chain id), or 9,999 residues per chain (wwPDB file-format documentation; the PDB archive itself is Berman et al. 2000 *Nucleic Acids Res* 28:235-242). Large assemblies (ribosomes, capsids, spliceosomes, most big cryo-EM structures) therefore exist ONLY as mmCIF, and converting a big mmCIF down to PDB renames multi-character chains and overflows serial numbers, silently corrupting any downstream tool that keys on chain id. mmCIF has been the wwPDB archive standard since 2014 and mandatory for crystallographic depositions since July 2019.

Two further traps compound this. First, Bio.PDB is PERMISSIVE by design: it reads malformed files, and it silently drops anisotropic B-factors (ANISOU), collapses each disordered atom to its highest-occupancy alternate, and never models most metadata (resolution, method, R-free, entity graph, assembly operators). Parse-success is not data integrity. Second, the deposited coordinates for an X-ray entry are usually the ASYMMETRIC UNIT, a crystallographic bookkeeping object that is frequently NOT the biologically functional oligomer -- so any interface, oligomeric-state, or buried-surface question must first obtain the biological assembly (Krissinel & Henrick 2007 *J Mol Biol* 372:774). "One chain in the file" is never evidence of a monomer.

The escape hatch for all three ceilings (assembly generation, very large structures, full mmCIF fidelity) is gemmi (Wojdyr 2022 *JOSS* 7:4200); Bio.PDB cannot apply the assembly operators itself. Prefer Bio.PDB for teaching, small structures, and hierarchy walks; reach for gemmi when the questions above appear.

## Decision: which format

| Format | Best when | Fails when | Hard limits |
|--------|-----------|------------|-------------|
| mmCIF / PDBx (`.cif`, `.cif.gz`) | Any modern default; large assemblies; full metadata; auth+label numbering; ANISOU/entities | A legacy tool only reads fixed-column PDB | None |
| Legacy PDB (`.pdb`, `.ent`) | Small structure feeding an old tool that demands PDB columns | Structure exceeds the format's limits (silently truncates/renames) | 99,999 atoms, ~62 chains, 9,999 resseq/chain, single-char chain id |
| BinaryCIF (`.bcif`, `.bcif.gz`) | Compact binary transport at bandwidth/scale; the current binary format | An ecosystem still expects the retired MMTF | None (lossless mmCIF encoding) |
| MMTF (`.mmtf`) | Nothing new -- RCSB stopped serving MMTF on 2 July 2024 | Any live download (the endpoint is decommissioned); treat as read-only-legacy | Retired upstream |

## Decision: Bio.PDB vs gemmi

| Task | Tool | Why |
|------|------|-----|
| Hierarchy walk, small X-ray/NMR structure, teaching | Bio.PDB | Readable SMCRA tree, pure Python, ubiquitous |
| Read metadata Bio.PDB drops (resolution, method, R-free, assembly ops) | Bio.PDB `MMCIF2Dict` | Raw category access without object-model loss |
| Generate the biological assembly from deposited coords | gemmi | Applies `_pdbx_struct_oper_list`; Bio.PDB has no operator-application code |
| Very large structure (>100k atoms), many structures, fast neighbor search | gemmi | C++ core scales; Bio.PDB's pure-Python tree is slow/memory-heavy |
| mmCIF round-trip without data loss (entities, label scheme, ANISOU) | gemmi | Full PDBx data model; writes hybrid-36 serials when >99,999 |

## Required Imports

```python
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, MMCIFIO, PDBList, Select
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.binary_cif import BinaryCIFParser
```

## Parse an mmCIF File (the modern default)

```python
from Bio.PDB import MMCIFParser

# auth_residues/auth_chains default to True: numbering matches the paper/UniProt.
parser = MMCIFParser(QUIET=True)
structure = parser.get_structure('4hhb', '4hhb.cif')

# label numbering is contiguous 1..N with no insertion codes -- a DIFFERENT scheme.
label_parser = MMCIFParser(QUIET=True, auth_residues=False, auth_chains=False)
label_structure = label_parser.get_structure('4hhb', '4hhb.cif')
```

Setting `auth_residues=False` renumbers to the mmCIF internal label scheme, so residue 100 in one parse is a different residue in the other. This is the single most common "my selection points at the wrong residue" bug. auth is what matches the literature and sequence databases; label is gap-free internal bookkeeping. Pick one scheme and stay in it.

## Parse a Legacy PDB File

```python
from Bio.PDB import PDBParser

# QUIET=True suppresses PDBConstructionWarning (discontinuous chains, missing occupancy).
parser = PDBParser(QUIET=True)
structure = parser.get_structure('1crn', '1crn.pdb')
```

## Parse a BinaryCIF File (compact binary, replaces MMTF)

```python
from Bio.PDB.binary_cif import BinaryCIFParser

# .get_structure(id, source); gz is handled transparently.
parser = BinaryCIFParser()
structure = parser.get_structure('1gbt', '1gbt.bcif.gz')
```

MMTF is intentionally absent here. RCSB retired it on 2 July 2024 and `MMTFParser.get_structure_from_url` targets a decommissioned service; BinaryCIF is the replacement.

## Read Metadata Bio.PDB Drops (MMCIF2Dict)

```python
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

# MMCIF2Dict returns category -> list[str]; index [0] and cast yourself.
meta = MMCIF2Dict('4hhb.cif')
resolution = meta.get('_refine.ls_d_res_high', ['NA'])[0]
method = meta.get('_exptl.method', ['NA'])[0]
r_free = meta.get('_refine.ls_R_factor_R_free', ['NA'])[0]
r_work = meta.get('_refine.ls_R_factor_R_work', ['NA'])[0]
```

The parser's thin `structure.header` omits resolution/R-free for many files; the dict reaches anything in the mmCIF, including assembly operators the object model never builds.

## Download from RCSB (PDBList)

```python
from Bio.PDB import PDBList

pdbl = PDBList()

# pdir=None writes into a two-char divided subdirectory tree (e.g. hh/4hhb.cif),
# NOT the current directory; pass pdir='.' to control the location.
path = pdbl.retrieve_pdb_file('4HHB', pdir='.', file_format='mmCif')

# file_format 'pdb' fetches legacy PDB only when the entry fits the format.
legacy_path = pdbl.retrieve_pdb_file('4HHB', pdir='.', file_format='pdb')
```

`file_format='mmCif'` (that exact casing) is the current default recommendation. `retrieve_pdb_file` has no `assembly_num` parameter in current Biopython -- download the assembly directly (below).

## Download the Biological Assembly (Bio.PDB cannot build it)

```python
import gzip, shutil, urllib.request

# The ASU is often not the functional oligomer; RCSB pre-applies the operators
# in the -assemblyN file, so downloading it is safer than regenerating.
pdb_id = '1abc'
url = f'https://files.rcsb.org/download/{pdb_id.upper()}-assembly1.cif.gz'
urllib.request.urlretrieve(url, f'{pdb_id}-assembly1.cif.gz')
with gzip.open(f'{pdb_id}-assembly1.cif.gz', 'rb') as fin, open(f'{pdb_id}-assembly1.cif', 'wb') as fout:
    shutil.copyfileobj(fin, fout)
```

Bio.PDB has no operator-application code, so if only the deposited ASU is on disk it cannot construct the assembly -- use the RCSB assembly file or gemmi's `transform_to_assembly`.

## Write a Structure

```python
from Bio.PDB import PDBParser, MMCIFIO

parser = PDBParser(QUIET=True)
structure = parser.get_structure('1crn', '1crn.pdb')

# Writing mmCIF preserves multi-char chains and >99,999 serials; PDB cannot.
io = MMCIFIO()
io.set_structure(structure)
io.save('1crn_out.cif')
```

## Write a Subset with the Select Class

```python
from Bio.PDB import PDBParser, PDBIO, Select

class ProteinChainSelect(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        # id[0] is the hetflag: ' ' standard, 'W' water, 'H_XXX' hetero.
        return residue.id[0] == ' '

parser = PDBParser(QUIET=True)
structure = parser.get_structure('1crn', '1crn.pdb')

io = PDBIO()
io.set_structure(structure)
io.save('chain_A_protein.pdb', ProteinChainSelect('A'))
```

Override any of `accept_model`, `accept_chain`, `accept_residue`, `accept_atom` to return truthy to keep. Writing a large mmCIF back out as PDB through PDBIO is where multi-character chains and serial overflow silently corrupt the output.

## Capture Warnings for an Unfamiliar File

```python
from Bio.PDB import PDBParser
import warnings

# QUIET=True is the reflex, but it hides 'chain is discontinuous' -- the warning
# that flags a numbering gap or a merge the parser should not have made.
parser = PDBParser(QUIET=False)
with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter('always')
    structure = parser.get_structure('unknown', 'unknown.pdb')
    for w in caught:
        print(w.message)
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Residue numbers do not match the paper / a UniProt mapping | Parsed with `auth_residues=False`, so numbering is the label scheme | Use the default `auth_residues=True`; only switch to label for gap-free internal indexing, never mix schemes |
| `resolution`/`R-free` is `None` from `structure.header` | Bio.PDB's header dict is thin and omits refinement metadata | Read `_refine.ls_d_res_high`, `_refine.ls_R_factor_R_free`, `_exptl.method` via `MMCIF2Dict` |
| `TypeError: retrieve_pdb_file() got an unexpected keyword 'assembly_num'` | Current Biopython `PDBList` has no `assembly_num` parameter | Download `...-assembly1.cif.gz` from files.rcsb.org directly (or use gemmi) |
| Downloaded file is not in the current directory | `pdir=None` writes a two-char divided subdirectory tree (`hh/4hhb.cif`) | Pass an explicit `pdir='.'` (or the target dir) to `retrieve_pdb_file` |
| MMTF download 404s / connection fails | RCSB retired MMTF on 2 July 2024; the endpoint is gone | Use BinaryCIF (`.bcif`) or mmCIF; treat MMTF files as read-only-legacy |
| Chains renamed and atom serials wrong after PDB output | A large mmCIF exceeded PDB's ~62-chain / 99,999-atom limits on write | Stay in mmCIF (`MMCIFIO`), or use gemmi's hybrid-36 writer |
| Analyzing a "monomer" that is really half a dimer | Computed on the deposited ASU, not the biological assembly | Fetch the `-assembly1` file (or generate with gemmi) before any interface/oligomer analysis |
| Download URL 404s for a newly deposited entry | The 4-char id space is being exhausted (~2028) and RCSB is phasing in extended 12-char ids (`pdb_00006uv8`) | Use the full extended id in the `files.rcsb.org` path; a hard-coded 4-char assumption breaks once extended ids arrive |
| Anisotropic B-factors (ANISOU) lost after a Bio.PDB round-trip | ANISOU is parsed but not reliably written back | Preserve the original file, or round-trip through gemmi when ANISOU matters |
| Distances/clashes look wrong at a partially disordered site | A disordered atom silently forwards to its highest-occupancy altloc | Enumerate altlocs with `atom.disordered_get_list()` and set an explicit altloc policy (see structure-navigation) |
| `KeyError` fetching a residue by integer, e.g. `chain[100]` | Residue id is the tuple `(hetflag, resseq, icode)`; insertion codes and hetero break the bare-int path | Key on the full 3-tuple, e.g. `chain[(' ', 100, ' ')]` |
| `BinaryCIFParser` import fails from `Bio.PDB` | It lives in the submodule `Bio.PDB.binary_cif`, not the top-level package | `from Bio.PDB.binary_cif import BinaryCIFParser` |
| Silent wrong results on a malformed file that "parsed fine" | Bio.PDB is permissive; parse-success is not data integrity | Parse with `QUIET=False` and inspect `PDBConstructionWarning`s for unfamiliar files |

## Related Skills

- structure-navigation - Walk the SMCRA tree, handle altlocs/insertion codes, extract observed vs SEQRES sequence
- structure-modification - Transform coordinates, strip waters/hetero safely, edit B-factors before writing
- geometric-analysis - Measure distances, angles, SASA, and superimpose once the correct assembly is loaded
- interface-analysis - Analyze the interfaces that only exist in the biological assembly, not the ASU
- structure-validation - Read resolution/R-free/clashscore to judge whether the loaded model is trustworthy
- structure-preparation - Add hydrogens/protonation and fill atoms on the loaded assembly before docking or MD
- alignment/structural-alignment - Superpose sequence-different structures that Bio.PDB's ordered correspondence cannot handle
- database-access/uniprot-access - Map structure residues back to a UniProt reference sequence

## References

- Berman HM, Westbrook J, Feng Z, Gilliland G, Bhat TN, Weissig H, Shindyalov IN, Bourne PE (2000). The Protein Data Bank. *Nucleic Acids Res* 28(1):235-242. DOI 10.1093/nar/28.1.235.
- Cock PJA, Antao T, Chang JT, Chapman BA, Cox CJ, Dalke A, Friedberg I, Hamelryck T, Kauff F, Wilczynski B, de Hoon MJL (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics* 25(11):1422-1423. DOI 10.1093/bioinformatics/btp163.
- Hamelryck T, Manderick B (2003). PDB file parser and structure class implemented in Python. *Bioinformatics* 19(17):2308-2310. DOI 10.1093/bioinformatics/btg332.
- Wojdyr M (2022). GEMMI: A library for structural biology. *Journal of Open Source Software* 7(73):4200. DOI 10.21105/joss.04200.
- Kunzmann P, Hamacher K (2018). Biotite: a unifying open source computational biology framework in Python. *BMC Bioinformatics* 19:346. DOI 10.1186/s12859-018-2367-z.
- Kim H, Mirdita M, Steinegger M (2023). Foldcomp: a library and format for compressing and indexing large protein structure sets. *Bioinformatics* 39(4):btad153. DOI 10.1093/bioinformatics/btad153.
- Krissinel E, Henrick K (2007). Inference of macromolecular assemblies from crystalline state. *J Mol Biol* 372(3):774-797. DOI 10.1016/j.jmb.2007.05.022.
- RCSB PDB (2024). Switch from MMTF to BinaryCIF: RCSB ceased serving MMTF on 2 July 2024. https://www.rcsb.org/news/65a1af31c76ca3abcc925d0c
- wwPDB. File formats and the PDB (legacy format frozen 2012; mmCIF archive standard 2014; 99,999-atom / 62-chain limits; large entries mmCIF-only). https://www.wwpdb.org/documentation/file-formats-and-the-pdb

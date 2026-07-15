---
name: bio-substructure-search
description: Searches molecular libraries for substructure matches using SMARTS patterns with explicit handling of recursive SMARTS, ring membership, aromaticity dialect, vector binding, atom map indices, and reactive/PAINS/REOS/Brenk filter catalogs. Use when filtering compounds by pharmacophore features, functional groups, scaffold matches, or screening for assay-interference / structural alerts.
tool_type: python
primary_tool: RDKit
---

## Version Compatibility

Reference examples tested with: RDKit 2024.09+. SMARTS dialect follows Daylight specification with RDKit extensions.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show rdkit` then `help(rdkit.Chem.MolFromSmarts)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Substructure Search

Search molecular collections for structural patterns using SMARTS. The choice of SMARTS dialect, atom/bond matching mode, and structural-alert catalog determines whether the search is correctly capturing the intended chemistry. PAINS (Baell & Holloway 2010) is the most-cited but most-misunderstood filter -- it identifies patterns of assay interference, not "bad molecules". Knowing when to apply each catalog and how to interpret hits is essential.

For SMARTS-based reactions (transforming matched substructures), see `chemoinformatics/reaction-enumeration`. For 3D pharmacophore matching, see `chemoinformatics/pharmacophore-modeling`.

## SMARTS Grammar Essentials

| Token | Meaning | Example |
|-------|---------|---------|
| `[#6]` | Atom by atomic number | `[#6]` carbon (any hybridization) |
| `c` | Lowercase = aromatic | `c1ccccc1` benzene aromatic |
| `C` | Uppercase = aliphatic only | `C(=O)O` carboxylic acid carbon |
| `[CX4]` | Atom + connection count X | `[CX4]` sp3 carbon (4 connections) |
| `[CX3]=O` | Carbonyl (CX3 = sp2 with 3 bonds) | matches ketone, aldehyde, ester C |
| `[#6;R]` | Atom in ring | `[#6;R]` ring carbon |
| `[#6;!R]` | Atom not in ring | `[#6;!R]` acyclic carbon |
| `[#6;r6]` | Atom in 6-membered ring | `[#6;r6]` six-ring carbon |
| `[a]` | Any aromatic atom | `[a]` |
| `[!#1]` | Anything except H | `[!#1]` heavy atom |
| `[N;H2]` | N with exactly 2 H; neighboring chemistry unconstrained | `[NH2]` also matches non-amine `NH2` environments unless context is added |
| `[N+]` | Positively charged N | `[N+](=O)[O-]` nitro |
| `[$(...)]` | Recursive SMARTS | `[$(c1ccccc1)]` aromatic 6-ring atom |
| `[c]([F,Cl,Br,I])` | OR within brackets | aryl halide |
| `~` | Any bond type | `c~c` any aromatic-aromatic bond |
| `@` | Ring bond | `c@c` requires the matched bond to be in a ring |
| `-` | Single bond explicit | `C-C` |
| `=` | Double bond | `C=O` |
| `:` | Aromatic bond explicit | |

## Common SMARTS Patterns

| Pattern | SMARTS | Notes |
|---------|--------|-------|
| Hydroxyl (alcohol + phenol) | `[OX2H]` | OX2H avoids matching O- in OH- |
| Phenol only | `[OX2H][c]` | OH attached to aromatic carbon |
| Aliphatic OH only | `[OX2H][CX4]` | OH attached to sp3 C |
| Carboxylic acid | `[CX3](=O)[OX2H1]` | C(=O)OH |
| Carboxylate | `[CX3](=O)[O-]` | C(=O)O- (deprotonated) |
| Ester | `[CX3](=O)[OX2][!H]` | C(=O)O-R |
| Amide | `[CX3](=[OX1])[NX3]` | C(=O)N-R |
| Primary amine attached to carbon (excluding common amide-like N) | `[NX3;H2;$(N-[#6]);!$(N-[C,S,P]=[O,S,N])]` | Carbon-substituted -NH2; extend the exclusions for a project-specific amine definition |
| Secondary amine attached to two carbons | `[NX3;H1;$(N(-[#6])-[#6]);!$(N-[C,S,P]=[O,S,N])]` | Carbon-substituted -NH- excluding common amide-like N |
| Neutral tertiary amine attached to three carbons | `[NX3;H0;+0;$(N(-[#6])(-[#6])-[#6]);!$(N-[C,S,P]=[O,S,N])]` | Carbon-substituted -NR2 excluding common amide-like N |
| Quaternary amine | `[NX4+]` | -NR4+ |
| Nitro | `[N+](=O)[O-]` | -NO2 |
| Nitrile | `[CX2]#[NX1]` | -C#N |
| Sulfonamide | `[SX4](=[OX1])(=[OX1])[NX3]` | -S(=O)(=O)N |
| Aryl halide | `[c][F,Cl,Br,I]` | halogen on aromatic |
| Aliphatic halide | `[CX4][F,Cl,Br,I]` | halogen on sp3 C |
| Hydrogen bond donor | Use a named feature definition such as RDKit `BaseFeatures.fdef` or `Lipinski.NumHDonors` | `[#7,#8;!H0]` is only a simplified N/O-H query and is not a universal HBD model |
| Hydrogen bond acceptor | Use a named feature definition such as RDKit `BaseFeatures.fdef` or `Lipinski.NumHAcceptors` | No short universal SMARTS correctly captures every accepted HBA chemistry model |
| Michael acceptor | `[CX3]=[CX3][CX3]=O` | enone, acrylamide warhead |
| Aldehyde | `[CX3H1](=O)` | -CHO |
| Ketone | `[CX3;H0](=[OX1])([#6])[#6]` | Carbonyl carbon has two carbon substituents and no hydrogen |

## Basic Substructure Match

**Goal:** Test whether a molecule contains a SMARTS pattern and enumerate the matching atom indices.

**Approach:** Parse the molecule with `MolFromSmiles` and the pattern with `MolFromSmarts`, gate with `HasSubstructMatch`, then call `GetSubstructMatches` and map each atom index back to the molecule for inspection.

```python
from rdkit import Chem

mol = Chem.MolFromSmiles('c1ccc(O)cc1CCO')
pattern = Chem.MolFromSmarts('[OX2H]')

if mol.HasSubstructMatch(pattern):
    matches = mol.GetSubstructMatches(pattern)
    for match in matches:
        atoms = [mol.GetAtomWithIdx(i).GetSymbol() for i in match]
```

`HasSubstructMatch` returns bool, `GetSubstructMatches` returns tuple of tuples of atom indices.

## Recursive SMARTS for Context-Aware Patterns

`[$(pattern)]` matches an atom that *also* matches the entire pattern starting from itself. Critical for context-aware matching.

```python
# Aromatic carbon attached to a carbonyl
pat = Chem.MolFromSmarts('[$(c[C](=O))]')

# Aniline-type N (aromatic carbon-N-H)
pat = Chem.MolFromSmarts('[$([NX3;H2][c])]')

# Neutral tertiary amine with three sp3-carbon neighbors
pat = Chem.MolFromSmarts('[$([NX3]([CX4])([CX4])[CX4])]')

# H-bond donor (per Lipinski, exclude quaternary)
hbd = Chem.MolFromSmarts('[#7,#8;!H0;!$([NX3+])]')

# For H-bond acceptors, use RDKit's maintained Lipinski/feature definitions
# instead of an ad hoc universal SMARTS.
from rdkit.Chem import Lipinski
n_acceptors = Lipinski.NumHAcceptors(mol)
```

## Structural-Alert Filter Catalogs

| Filter | Origin | Patterns | Use case | Failure mode |
|--------|--------|----------|----------|--------------|
| PAINS_A | Baell & Holloway 2010 | 16 | Most populated source-data patterns (>=150 analogues per pattern) | Many false positives in primary screens; legitimate medicines flagged |
| PAINS_B | Baell & Holloway 2010 | 55 | Intermediate source-data population (15-149 analogues per pattern) | Similar |
| PAINS_C | Baell & Holloway 2010 | 409 | Least populated source-data patterns (1-14 analogues per pattern) | Most permissive |
| BRENK | Brenk 2008 (DDS unsuitable) | 105 | Reactive / toxicity / undesirable | Useful for fragment / virtual library |
| NIH | NIH MLSMR | 180 in RDKit 2024.09 | Reactive groups, unstable | Legacy filter; verify count after toolkit upgrades |
| ZINC | ZINC clean-leads | 50 in RDKit 2024.09 | Drug-like cleanup | Verify definitions after toolkit upgrades |
| Glaxo / Eli Lilly | Vendor lists | varies | Internal "ugly" filters | Often unpublished |
| REOS | Walters & Murcko 2002 | property + structural | Drug-likeness combined filter | Hand-curated thresholds |

The PAINS A/B/C families encode pattern population in the original screening dataset, not increasing or decreasing external evidence strength.

## When to Apply Each Filter

| Scenario | Catalog | Reason |
|----------|---------|--------|
| Hit validation from biochemical screen | PAINS_A | Identify assay-interference candidates |
| Library prep for HTS | PAINS_A + Brenk + ZINC | Remove clearly bad |
| Fragment library design | Brenk + ZINC | Remove reactive; PAINS less critical at fragments |
| Lead optimization | None mandatory | Filters can exclude valid leads |
| Natural product analog | None | Filters trained on synthetic chemistry |
| Covalent inhibitor design | Skip warhead filter | Warheads ARE the design |

**Critical:** Capuzzi et al. (2017) found PAINS alerts in 87 FDA-approved small-molecule drugs. PAINS is a *flag for assay validation*, not a *killing filter*.

## PAINS Filter

**Goal:** Split a molecule list into PAINS-flagged and PAINS-clean sets using one or more PAINS catalog tiers.

**Approach:** Configure `FilterCatalogParams` with the requested catalog enums, build a `FilterCatalog` once, and for each molecule use `GetFirstMatch` to either bucket it as clean or record the matching pattern description.

```python
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

def pains_filter(mols, catalogs=('PAINS_A',)):
    params = FilterCatalogParams()
    for cat in catalogs:
        params.AddCatalog(getattr(FilterCatalogParams.FilterCatalogs, cat))
    catalog = FilterCatalog(params)

    flagged = []
    clean = []
    for mol in mols:
        if mol is None:
            continue
        entry = catalog.GetFirstMatch(mol)
        if entry is None:
            clean.append(mol)
        else:
            flagged.append((mol, entry.GetDescription()))
    return clean, flagged
```

Available catalog names: `PAINS_A`, `PAINS_B`, `PAINS_C`, `PAINS` (all), `BRENK`, `NIH`, `ZINC`, `ALL`.

## Reaction-Reactive Group Filter (custom)

For HTS triage, filter electrophilic warheads (acrylamide, chloroacetamide, etc.) unless designing covalent inhibitors.

**Goal:** Flag molecules containing electrophilic warheads or other reactive functional groups that would interfere with biochemical HTS.

**Approach:** Maintain a named SMARTS dictionary of reactive groups (acid halides, epoxides, Michael acceptors, etc.), then per molecule scan each pattern with `HasSubstructMatch` and return the first matching warhead name.

```python
REACTIVE_SMARTS = {
    'acid_anhydride': '[CX3](=O)O[CX3](=O)',
    'acid_halide': '[CX3](=O)[F,Cl,Br,I]',
    'alpha_halo_carbonyl': '[CX3](=O)C([F,Cl,Br,I])',
    'aldehyde_reactive': '[CX3H1](=O)[#6;X4]',  # aliphatic aldehydes
    'epoxide': 'C1OC1',
    'aziridine': 'C1NC1',
    'isocyanate': '[NX2]=C=[OX1]',
    'isothiocyanate': '[NX2]=C=[SX1]',
    'beta_lactam': 'C1(=O)NCC1',
    'sulfonyl_halide': '[SX4](=O)(=O)[F,Cl,Br,I]',
    'Michael_acceptor': '[CX3]=[CX3][CX3]=O',
    'vinyl_sulfone': '[SX4](=O)(=O)C=C',
}

def reactive_filter(mol, exclude_warheads=True):
    if not exclude_warheads:
        return False, None
    for name, smarts in REACTIVE_SMARTS.items():
        if mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
            return True, name
    return False, None
```

For covalent-inhibitor design, see `chemoinformatics/covalent-design`; these warheads are the desired chemistry, not noise to filter.

## Library Filtering with Multiple Patterns

**Goal:** Reduce a molecule library to those that match all required SMARTS patterns and none of the excluded ones.

**Approach:** Start from the full molecule list, iteratively intersect with each `include` SMARTS using `HasSubstructMatch`, then subtract any molecule matching an `exclude` SMARTS.

```python
def filter_library(mols, include=None, exclude=None):
    keep = list(mols)
    if include:
        for s in include:
            p = Chem.MolFromSmarts(s)
            keep = [m for m in keep if m and m.HasSubstructMatch(p)]
    if exclude:
        for s in exclude:
            p = Chem.MolFromSmarts(s)
            keep = [m for m in keep if m and not m.HasSubstructMatch(p)]
    return keep
```

## Atom Map Indices in SMARTS

Atom maps `[C:1]` track atoms through transformations. Used in reactions (`reaction-enumeration` skill) but also for substructure-based extraction:

```python
# Find amide N with attached aryl
pat = Chem.MolFromSmarts('[CX3:1](=O)[NX3:2][c:3]')
match = mol.GetSubstructMatch(pat)

amide_C, amide_N, aryl_C = match
```

## Per-Tool Failure Modes

### PAINS -- false positive on natural product

**Trigger:** Library contains natural products, polyphenols, flavonoids, quinones.

**Mechanism:** PAINS_A patterns target rhodanines, curcumins, polyhydroxylated polyphenols -- legitimate scaffolds in natural-product chemistry.

**Symptom:** Library hits flagged as PAINS but trace back to validated natural products with confirmed activity.

**Fix:** Use PAINS as a *flag* not a *delete*. Cross-check flagged compounds for orthogonal-assay confirmation (label-free e.g. SPR, ITC).

### Aromaticity dialect mismatch

**Trigger:** SMARTS pattern with `c` (aromatic) for a heteroatom-rich ring; molecule parsed with different aromaticity model.

**Mechanism:** RDKit, OpenEye, ChemAxon differ on whether furan, thiazole, tropone, etc. are aromatic.

**Symptom:** Same pattern matches in one toolkit, not in another.

**Fix:** Re-canonicalize molecules within RDKit before applying SMARTS. Or use `[#6]:[#6]` instead of `c:c` (explicit element + bond type).

### Tautomer-sensitive pattern miss

**Trigger:** SMARTS targets keto form `C(=O)` but molecule is enol `C(O)=C`.

**Mechanism:** Default canonical form differs by toolkit + standardization choice.

**Symptom:** Known matching molecule reports no match.

**Fix:** Use tautomer-aware match: enumerate tautomers and OR-match. Or canonicalize first via `chemoinformatics/molecular-standardization`. Or expand pattern with `[$(C(=O)),$(C(O)=C)]`.

### Stereochemistry ignored

**Trigger:** SMARTS without `/\@` stereo markers applied to mol with explicit stereo.

**Mechanism:** SMARTS matching is stereo-agnostic by default.

**Symptom:** Wrong stereoisomer is matched as well as right one.

**Fix:** `mol.GetSubstructMatches(pattern, useChirality=True)` to require chirality match.

### Ring closure / fused-ring specificity

**Trigger:** A query must distinguish an isolated benzene ring from a six-membered aromatic ring embedded in a fused system.

**Mechanism:** `c1ccccc1` matches six-membered aromatic cycles and therefore does match benzene cycles within naphthalene. Extra ring-membership or fusion constraints are required to exclude fused systems.

**Symptom:** A nominal "benzene" query returns fused polyaromatics that the project intended to exclude.

**Fix:** Keep `c1ccccc1` when any aromatic six-cycle is desired. When an isolated ring is required, add explicit ring-degree/fusion constraints and test the query against benzene, naphthalene, indole, and representative substituted controls.

### Recursive SMARTS performance

**Trigger:** Deeply nested recursive SMARTS over a large library.

**Mechanism:** Each `[$()]` re-evaluates the inner pattern for every candidate atom.

**Symptom:** Search 10x-100x slower than expected.

**Fix:** Flatten recursion where possible; pre-filter with simpler pattern, then re-test with the recursive one.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| `Chem.MolFromSmarts` returns None | Invalid SMARTS grammar | Validate with `Chem.MolFromSmarts(smi, mergeHs=False)`; check parens, brackets |
| `[OH]` gives unexpected hydroxyl matches | Query does not state the intended valence/connectivity model | Use `[OX2H]` for neutral alcohol/phenol oxygen or a more specific context-aware pattern |
| Pattern matches but library is "empty" | Mol failed sanitize | Try `Chem.SDMolSupplier(sanitize=False)` then catch errors |
| Multiple matches per molecule | Single-match query expected | `GetSubstructMatch` returns first; `GetSubstructMatches` returns all |
| Match indices but no fragment | Match returns atom indices in pattern order | Map to original mol via `mol.GetAtomWithIdx(i)` |
| PAINS catalog initialization slow | Loading 1000+ patterns on every call | Build catalog once, reuse for batch |
| Stereo SMARTS not matching | `useChirality=False` (default) | `mol.GetSubstructMatches(p, useChirality=True)` |

## References

- Baell & Holloway, *J. Med. Chem.* 53:2719-2740 (2010) -- original PAINS filter and tier evidence. https://doi.org/10.1021/jm901137j
- Capuzzi et al., *J. Chem. Inf. Model.* 57:417-427 (2017) -- PAINS reality check (FDA drug overlap). https://doi.org/10.1021/acs.jcim.6b00465
- Brenk et al., *ChemMedChem* 3:435-444 (2008) -- structural alerts (BRENK filter). https://doi.org/10.1002/cmdc.200700139
- Walters & Murcko, *Adv. Drug Deliv. Rev.* 54:255-271 (2002) -- drug-likeness filtering, including REOS. https://doi.org/10.1016/S0169-409X(02)00003-0
- Bruns & Watson, *J. Med. Chem.* 55:9763-9772 (2012) -- Eli Lilly medchem rules. https://doi.org/10.1021/jm301008n
- Daylight Chemical Information Systems, SMARTS theory documentation -- complete grammar reference. https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html

## Related Skills

- chemoinformatics/molecular-io - Parse molecules before searching
- chemoinformatics/molecular-standardization - Canonicalize tautomers before SMARTS
- chemoinformatics/similarity-searching - Fingerprint-based fuzzy matching
- chemoinformatics/scaffold-analysis - Scaffold-based pattern derivation
- chemoinformatics/reaction-enumeration - SMARTS for chemical transformations
- chemoinformatics/admet-prediction - PAINS as ADMET filter
- chemoinformatics/covalent-design - Warhead chemistry

---
name: bio-covalent-design
description: Designs covalent inhibitors and warheads targeting cysteine, lysine, serine, threonine, tyrosine, and aspartate residues, with explicit handling of warhead reactivity (acrylamide, chloroacetamide, vinyl sulfone, sulfonyl fluoride, fluorosulfate, aldehyde, boronate, nitrile), reversibility (kinact/Ki, t_residence), glutathione (GSH) stability, intrinsic reactivity assays, and covalent docking (DOCKovalent, GOLD, HCovDock). Use when designing covalent inhibitors for targeted covalent inhibition (TCI), KRAS G12C-style approaches, or rationalizing covalent SAR.
tool_type: python
primary_tool: RDKit
---

## Version Compatibility

Reference examples tested with: RDKit 2024.09+, OpenEye / AutoDock Vina 1.2+ (for covalent extensions), GOLD (commercial), DOCKovalent (web service), HCovDock 1.0+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show rdkit` then `help(rdkit.Chem)` to check signatures
- CLI: check version output of each docking tool

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Covalent Inhibitor Design

Design molecules that form covalent bonds with target protein residues. Clinically validated targeted covalent inhibitors include KRAS G12C inhibitors (sotorasib, adagrasib), BTK inhibitors (ibrutinib), and EGFR inhibitors (osimertinib). Covalent design requires balancing **intrinsic reactivity** (must form bond) vs **selectivity** (only the intended residue), **reversibility** (irreversible vs reversible covalent), and **drug-likeness** (warheads can hurt PK).

For warhead substructure filtering (in non-covalent contexts), see `chemoinformatics/substructure-search`. For non-covalent docking, see `chemoinformatics/virtual-screening`. For pose validation, see `chemoinformatics/pose-validation`.

## Reactive Residue Taxonomy

| Residue | Nucleophile | Example compatible warheads | Design note |
|---------|-------------|-----------------------------|-------------|
| Cysteine | Thiol / thiolate | Acrylamide, haloacetamide, nitrile | Commonly targeted; local pKa and geometry control reactivity |
| Lysine | Amine | Sulfonyl fluoride, aldehyde | Aldehydes can form reversible imines with amines |
| Serine | Alcohol / alkoxide | β-lactam, boronate | Often requires catalytic activation |
| Threonine | Alcohol / alkoxide | Boronate | Context-dependent and less commonly targeted |
| Tyrosine | Phenol / phenolate | Sulfonyl fluoride, fluorosulfate | Local environment strongly affects reaction |
| Aspartate/Glutamate | Carboxylate | Residue-specific electrophiles require experimental validation | Do not infer aldehyde Schiff-base formation with carboxylates |

Cysteine is frequently targeted because its thiol/thiolate can be nucleophilic and its local environment can support selective proximity-driven reaction. GSH and off-target cysteines are competing thiols, not intrinsically distinguishable from the target by the warhead alone; the complete ligand's recognition, exposure, and intrinsic reactivity determine selectivity.

## Warhead Chemistry

| Warhead | SMARTS pattern | Reactivity | Reversibility | Cys-selective |
|---------|----------------|------------|---------------|----------------|
| Acrylamide | `[CX3](=[OX1])([NX3])[CX3]=[CX3]` | Moderate (Michael acceptor) | Usually irreversible | Often Cys-directed |
| Chloroacetamide | `[CX3](=[OX1])([NX3])[CH2]Cl` | High (SN2) | Irreversible | Often Cys-directed |
| α-haloketone | `[CX3](=O)C[F,Cl,Br]` | Very high | Irreversible | Yes (but reactive) |
| Vinyl sulfone | `S(=O)(=O)C=C` | Moderate (Michael) | Irreversible | Yes |
| Sulfonyl fluoride | `S(=O)(=O)F` | Moderate | Irreversible | Lys/Tyr/Ser |
| Fluorosulfate (SuFEx) | `OS(=O)(=O)F` | Moderate | Irreversible | Tyr/Lys |
| Aldehyde | `[CX3H1](=O)` | Variable | Often reversible (covalent equilibrium) | Context-dependent Cys/Lys/Ser chemistry |
| Boronate (B-OH or B(OH)2) | `B(O)O` | Moderate | Reversible | Ser/Thr |
| Nitrile | `C#N` | Low | Reversible (Cys-S adduct) | Cys |
| Epoxide | `C1OC1` | High | Irreversible | Cys/Lys/Asp |
| α,β-unsaturated ketone | `[CX3](=O)C=C` | Moderate (Michael) | Irreversible | Cys |
| Isothiocyanate | `N=C=S` | High | Irreversible | Cys/Lys |
| Maleimide | `O=C1N(C(=O)C=C1)` | Often high | Commonly irreversible under assay conditions | Often Cys-directed |
| Cysteine-selective heterocycle | various | Moderate | Variable | Yes (designed) |

**Practical hierarchy:** Acrylamides are common attenuated electrophiles in cysteine-directed TCIs, including KRAS G12C, EGFR, and BTK programs. Haloacetamides are generally more intrinsically reactive, but actual selectivity must be measured for the complete molecule and target context.

## Decision Tree by Scenario

| Goal | Warhead choice | Reactivity tier |
|------|----------------|-----------------|
| Cysteine TCI program | Acrylamide is one common starting class | Measure complete-compound target and intrinsic reactivity |
| Cysteine probe program | Haloacetamides are one common class | Higher intrinsic reactivity can aid labeling but requires selectivity profiling |
| Lysine TCI (uncommon) | Sulfonyl fluoride | Moderate |
| Tyrosine TCI | Fluorosulfate (SuFEx) | Moderate |
| Reversible covalent | Warhead with demonstrated reversible adduct chemistry, such as selected cyanoacrylamides, aldehydes, nitriles, or boronates | Confirm reversibility experimentally |
| Activity-based protein profiling (ABPP) | Iodoacetamide / chloroacetamide | Very high |
| Boronic acid inhibitor (proteasome) | Boronate | Reversible |
| Aldehyde inhibitor (calpain) | Aldehyde | Reversible covalent |

## Kinetics: kinact / Ki

Covalent inhibition kinetics:
- **Ki**: reversible binding affinity (initial, like non-covalent IC50)
- **kinact**: rate of covalent bond formation (sec^-1)
- **kinact/Ki**: second-order rate constant, "covalent efficiency" (M^-1 s^-1)

For irreversible two-step inhibition that follows the corresponding kinetic model, report fitted `kinact`, `Ki`, and `kinact/Ki` rather than only a time-dependent IC50. Two compounds with the same IC50 can have different kinetic components:
- Low Ki, low kinact: tight binding, slow covalent bond
- High Ki, high kinact: loose binding, fast covalent bond

Because `kinact/Ki` depends on the target, construct, assay conditions, and kinetic model, compare values within a matched assay series and alongside exposure, intrinsic reactivity, target engagement, and selectivity. Reversible-covalent systems may require different mechanistic models and residence-time or washout measurements.

## Intrinsic Reactivity Assays

Before committing to a warhead, measure intrinsic reactivity (off-target risk):

```python
# Generic GSH stability assay readout - measure half-life of warhead with 10 mM GSH
# kinact_GSH from time-course of warhead disappearance
```

Do not assign a universal GSH half-life from the warhead name alone. Substitution, electronics, ionization, solubility, and assay conditions can change the observed rate. Report the GSH concentration, buffer, temperature, analytical method, and fitted half-life or second-order rate constant for the complete compound.

## Covalent Docking Tools

| Tool | Approach | Strength | Fails when |
|------|----------|----------|------------|
| DOCKovalent (London et al 2014 Nat Chem Biol 10:1066) | Constraint-based DOCK | Free, well-validated | Browser-based; small library |
| GOLD covalent (CCDC) | GOLD with covalent constraint | Commercial; selectivity | License cost |
| AutoDock 4 covalent | AD4 with covalent bond | Open source | Slower than Vina |
| CovDock (Schrödinger) | Glide-based + covalent | Commercial two-stage covalent docking workflow | License cost |
| MOE covalent | Triposite Discovery | Commercial | License cost |
| HCovDock (Wu Q, Huang S-Y 2023 Briefings Bioinform 24:bbac559) | Hierarchical fragment + covalent | Open; supports many residues | Newer, less validated |
| ICM-Pro covalent | Active site grid + covalent | Commercial; metal centers | License cost |

For open-source covalent docking, **HCovDock** (2023) is the modern alternative; **DOCKovalent** is the longstanding standard.

## Example: KRAS G12C Inhibitor Design Workflow

**Goal:** Decorate a co-crystal scaffold with a cysteine-targeting warhead and rank candidates by covalent efficiency.

**Approach:** Load scaffold SMILES, enumerate acrylamide-bearing analogs, filter by reactivity selectivity, dock under covalent constraint, and rank by kinact/Ki surrogates.

```python
from rdkit import Chem

# Step 1: scaffold from co-crystal (4LRW or AMG510)
scaffold_smi = 'c1ccc(C(=O)NCC)cc1'  # generic valid scaffold for code illustration
scaffold = Chem.MolFromSmiles(scaffold_smi)

# Step 2: enumerate analogs with acrylamide warhead
def add_acrylamide(scaffold, attachment_atom_idx):
    """Project hook: attach a mapped acrylamide with an audited reaction."""
    raise NotImplementedError(
        'Provide a project-specific mapped reaction and validate atom mapping, '
        'valence, regioisomer identity, and product sanitization.'
    )

# Step 3: filter for reactive group selectivity
# Step 4: dock with DOCKovalent / GOLD covalent / HCovDock
# Step 5: rank by kinact/Ki surrogate (compute reactive Michael acceptor reactivity)
```

## Reactivity Surrogates (computed without experiment)

For ranking warheads without wet-lab data:

| Descriptor | Use case |
|------------|----------|
| LUMO energy (DFT) | Michael acceptor reactivity (lower LUMO = more reactive) |
| Electrophile partial charge | SN2 reactivity |
| RDKit `rdMolDescriptors.CalcLabuteASA` | Steric accessibility |
| Experimentally supported binding pose or validated docking model | Geometric fit to reactive residue |

**Goal:** Record alpha-carbon substitution as a structural feature for an acrylamide series.

**Approach:** Parse the SMILES, locate the acrylamide substructure, and count neighbors on the alpha carbon outside the matched warhead. This count is not a LUMO estimate or a stand-alone reactivity prediction; substituent electronics and the rest of the molecule must be considered, and reactivity should be measured.

```python
def acrylamide_alpha_substitution_count(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    acryl_pat = Chem.MolFromSmarts(
        '[CX3:1](=[OX1:2])([NX3:3])[CX3:4]=[CX3:5]'
    )
    matches = mol.GetSubstructMatches(acryl_pat, uniquify=True)
    if not matches:
        return None
    alpha_query_idx = next(
        atom.GetIdx() for atom in acryl_pat.GetAtoms()
        if atom.GetAtomMapNum() == 4
    )
    alpha_c = mol.GetAtomWithIdx(matches[0][alpha_query_idx])
    n_subs = len([n for n in alpha_c.GetNeighbors() if n.GetIdx() not in matches[0]])
    return n_subs
```

For a reactivity model, use experimentally measured rates or a validated quantum-chemical workflow; a single frontier-orbital energy is not sufficient on its own.

## Per-Tool Failure Modes

### Wrong warhead for residue

**Trigger:** A warhead/residue pairing is assumed from a broad class label.

**Mechanism:** Reaction depends on the residue microenvironment, electrophile, binding pose, and catalytic assistance; a class label does not establish residue selectivity.

**Symptom:** No covalent adduct observed despite docking pose.

**Fix:** Use literature-supported residue/warhead chemistry as a hypothesis, then verify site-specific adduct formation and competing reactivity experimentally.

### Excessive reactivity (off-target)

**Trigger:** Chloroacetamide in drug-candidate context.

**Mechanism:** Excess intrinsic electrophile reactivity can increase reaction with GSH and off-target nucleophiles.

**Symptom:** Toxicity in cell-based assays; non-specific binding signal.

**Fix:** Test a less intrinsically reactive electrophile and measure its GSH and target-reaction kinetics; alpha substitution can tune behavior but does not guarantee selectivity or stability.

### Geometric mismatch

**Trigger:** The ligand's reactive atom is poorly positioned relative to the target nucleophile.

**Mechanism:** Covalent reaction requires warhead-specific distance and approach geometry between the electrophilic atom and the nucleophilic atom (Cys Sγ for cysteine).

**Symptom:** No covalent labeling in mass spec despite predicted docking.

**Fix:** Identify the reaction atoms, inspect the pre-reaction Sγ-to-electrophile distance and reaction-specific angles, and use a docking protocol parameterized for that reaction. Do not substitute Cβ distance for the reacting sulfur.

### Reversibility unintended

**Trigger:** Designed irreversible TCI but warhead is reversible.

**Mechanism:** Reversibility depends on the complete electrophile, adduct chemistry, protein environment, and assay timescale; class-level labels are only hypotheses.

**Symptom:** Activity wanes after substrate washout in cellular assays.

**Fix:** Select chemistry with demonstrated behavior in the intended context and verify reversibility by dilution, washout, intact-protein MS, or another suitable experiment.

### kinact/Ki conflation

**Trigger:** Optimizing for IC50 instead of kinact/Ki.

**Mechanism:** Compounds with same IC50 differ in covalent efficiency.

**Symptom:** Compounds with similar endpoint IC50 values show different time-dependent target engagement or pharmacodynamic duration.

**Fix:** Fit a mechanistically appropriate kinetic model. Use `kinact/Ki` for qualifying irreversible two-step systems; use equilibrium, residence-time, or washout measurements where appropriate for reversible covalent systems.

### DOCKovalent over-prediction

**Trigger:** Default DOCKovalent run.

**Mechanism:** Covalent constraint forces docking; many ligands "succeed" but are unrealistic.

**Symptom:** Many compounds pass docking; few label in vitro.

**Fix:** Review measured/validated reactivity, reaction-atom geometry (Cys Sγ for cysteine), non-covalent recognition, strain, and site-specific experimental labeling.

## Reconciliation: Irreversible vs Reversible Covalent

| Aspect | Irreversible | Reversible covalent |
|--------|--------------|---------------------|
| Examples | KRAS G12C (acrylamide), BTK (ibrutinib) | Boronate (bortezomib), aldehyde (calpain inhibitors) |
| Toxicity profile | Off-target Cys labeling potential | Off-target equilibrium |
| Resistance mechanism | Mutation of reactive Cys | Mutation reduces affinity |
| Design decision | Consider duration of target engagement, safety, exposure, and resistance | Consider equilibrium, residence time, and recovery after washout |

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Warhead not matching SMARTS | Different stereochemistry or charged | Use canonicalized + neutral SMARTS |
| DOCKovalent rejects ligand | No suitable Cys in pocket | Re-check residue accessibility |
| GSH adduct dominates | Warhead too reactive | Use less reactive warhead; or alpha-substitute |
| Off-target labeling in cells | Promiscuous warhead | Iterate warhead reactivity vs selectivity |
| Docking pose but no labeling | Geometric mismatch | Distance check; rotamer search |
| Intended irreversible inhibitor shows recovery after washout | Adduct chemistry is reversible or covalent reaction is incomplete | Re-evaluate the mechanism and fit the appropriate kinetic model |
| HCovDock fails on PROTAC | Tool optimized for monomer covalent | Use specialized tools for bivalent |

## References

- Lonsdale & Ward, *Chem. Soc. Rev.* 47:3816-3830 (2018) -- irreversible-inhibitor discovery, optimization, and kinetics (DOI 10.1039/C7CS00720C).
- Singh J, Petter RC, Baillie TA, Whitty A. *Nat. Rev. Drug Discov.* 10:307-317 (2011) -- TCI design principles (DOI 10.1038/nrd3410).
- London N et al., *Nat. Chem. Biol.* 10:1066-1072 (2014) -- DOCKovalent (DOI 10.1038/nchembio.1666).
- Wu Q et al., *Brief. Bioinform.* 24:bbac559 (2023) -- HCovDock (DOI 10.1093/bib/bbac559).
- Yu W, Weber DJ, MacKerell AD Jr. *J. Chem. Theory Comput.* 19:3007-3021 (2023) -- SILCS-Covalent and Cys-sulfur/reactive-atom geometry (DOI 10.1021/acs.jctc.3c00232).
- Backus et al., *Nature* 534:570-574 (2016) -- proteome-wide covalent ligand discovery (DOI 10.1038/nature18002).
- Pettinger et al., *Angew. Chem. Int. Ed.* 56:15200-15209 (2017) -- lysine-targeting covalent inhibitors (DOI 10.1002/anie.201707630).
- Ostrem et al., *Nature* 503:548-551 (2013) -- KRAS G12C disulfide-tethered fragments (DOI 10.1038/nature12796).

## Related Skills

- chemoinformatics/molecular-io - Parse warhead SMILES
- chemoinformatics/substructure-search - Warhead SMARTS detection
- chemoinformatics/virtual-screening - Pre-dock candidate non-covalent fit
- chemoinformatics/pose-validation - Validate covalent docking
- chemoinformatics/conformer-generation - Warhead conformer ensembles
- chemoinformatics/admet-prediction - ADMET of covalent leads
- chemoinformatics/molecular-descriptors - Reactivity surrogate descriptors

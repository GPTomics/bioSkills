---
name: bio-covalent-design
description: Designs covalent inhibitors and warheads targeting cysteine (most common, 98% of covalent drugs), lysine, serine, threonine, tyrosine, and aspartate residues, with explicit handling of warhead reactivity (acrylamide, chloroacetamide, vinyl sulfone, sulfonyl fluoride, fluorosulfate, aldehyde, boronate, nitrile), reversibility (kinact/Ki, t_residence), glutathione (GSH) stability, intrinsic reactivity assays, and covalent docking (DOCKovalent, GOLD, HCovDock). Use when designing covalent inhibitors for targeted covalent inhibition (TCI), KRAS G12C-style approaches, or rationalizing covalent SAR.
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

Design molecules that form covalent bonds with target protein residues. The "covalent revolution" (Lonsdale & Ward 2018) made TCIs (Targeted Covalent Inhibitors) clinically validated: KRAS G12C inhibitors (sotorasib, adagrasib), BTK inhibitors (ibrutinib), and EGFR inhibitors (osimertinib) are recent successes. Postdoc-grade covalent design requires balancing **intrinsic reactivity** (must form bond) vs **selectivity** (only the intended residue), **reversibility** (irreversible vs reversible covalent), and **drug-likeness** (warheads can hurt PK).

For warhead substructure filtering (in non-covalent contexts), see `chemoinformatics/substructure-search`. For non-covalent docking, see `chemoinformatics/virtual-screening`. For pose validation, see `chemoinformatics/pose-validation`.

## Reactive Residue Taxonomy

| Residue | % of covalent drugs | Reactivity | Notes |
|---------|---------------------|------------|-------|
| Cysteine | ~98% | High (nucleophile thiol) | Most accessible; preferred |
| Lysine | ~1% | Moderate (amine) | Less reactive; selective for sulfonyl fluoride |
| Serine | <1% | Low (alcohol, requires activation) | β-lactam, boronate |
| Threonine | very rare | Low | Boronate, aldehyde |
| Tyrosine | very rare | Moderate (phenol) | Sulfonyl fluoride, fluorosulfate |
| Aspartate/Glutamate | very rare | Low (carboxylate) | Aldehyde Schiff base |

Cysteine is the dominant target because:
- Soft nucleophile (matches soft electrophiles)
- Low background reactivity (rare in proteins, ~1.7%)
- Distinguishable from common nucleophiles (GSH, off-target Cys)

## Warhead Chemistry

| Warhead | SMARTS pattern | Reactivity | Reversibility | Cys-selective |
|---------|----------------|------------|---------------|----------------|
| Acrylamide | `C(=O)C=C` | Moderate (Michael acceptor) | Irreversible | Yes |
| Chloroacetamide | `C(=O)CCl` | High (SN2) | Irreversible | Yes |
| α-haloketone | `[CX3](=O)C[F,Cl,Br]` | Very high | Irreversible | Yes (but reactive) |
| Vinyl sulfone | `S(=O)(=O)C=C` | Moderate (Michael) | Irreversible | Yes |
| Sulfonyl fluoride | `S(=O)(=O)F` | Moderate | Irreversible | Lys/Tyr/Ser |
| Fluorosulfate (SuFEx) | `OS(=O)(=O)F` | Moderate | Irreversible | Tyr/Lys |
| Aldehyde | `C(=O)[H]` | Variable | Reversible (covalent equilibrium) | Cys/Lys/Ser |
| Boronate (B-OH or B(OH)2) | `B(O)O` | Moderate | Reversible | Ser/Thr |
| Nitrile | `C#N` | Low | Reversible (Cys-S adduct) | Cys |
| Epoxide | `C1OC1` | High | Irreversible | Cys/Lys/Asp |
| α,β-unsaturated ketone | `[CX3](=O)C=C` | Moderate (Michael) | Irreversible | Cys |
| Isothiocyanate | `N=C=S` | High | Irreversible | Cys/Lys |
| Maleimide | `C(=O)N(C(=O))C=C` | Very high | Irreversible | Cys |
| Cysteine-selective heterocycle | various | Moderate | Variable | Yes (designed) |

**Practical hierarchy:** Acrylamide is the modern default for cysteine-selective TCIs (KRAS G12C, EGFR, BTK). Chloroacetamide is more reactive (faster) but less selective.

## Decision Tree by Scenario

| Goal | Warhead choice | Reactivity tier |
|------|----------------|-----------------|
| Cysteine TCI, drug candidate | Acrylamide | Moderate (~kinact/Ki ~10^3-10^5 M^-1 s^-1) |
| Cysteine probe (chemical biology) | Chloroacetamide | High (kinact/Ki ~10^4-10^6 M^-1 s^-1) |
| Lysine TCI (uncommon) | Sulfonyl fluoride | Moderate |
| Tyrosine TCI | Fluorosulfate (SuFEx) | Moderate |
| Reversible covalent (KRAS G12C-like) | Acrylamide with α-substitution | Moderate reversibility |
| Activity-based protein profiling (ABPP) | Iodoacetamide / chloroacetamide | Very high |
| Boronic acid inhibitor (proteasome) | Boronate | Reversible |
| Aldehyde inhibitor (calpain) | Aldehyde | Reversible covalent |

## Kinetics: kinact / Ki

Covalent inhibition kinetics:
- **Ki**: reversible binding affinity (initial, like non-covalent IC50)
- **kinact**: rate of covalent bond formation (sec^-1)
- **kinact/Ki**: second-order rate constant, "covalent efficiency" (M^-1 s^-1)

Modern best practice: report kinact/Ki, not just IC50. Two compounds with same IC50 can have very different kinact/Ki:
- Low Ki, low kinact: tight binding, slow covalent bond
- High Ki, high kinact: loose binding, fast covalent bond

| kinact/Ki range (M^-1 s^-1) | Interpretation | Reference |
|-----------------------------|----------------|-----------|
| > 10^5 | Highly efficient covalent inhibitor | Chloroacetamide probes, fragment-warhead TCIs |
| 10^3 - 10^5 | Standard for TCI; clinical candidate | Sotorasib (AMG510) KRAS G12C ~2x10^4 (Hallin 2020); adagrasib ~5x10^3 |
| 10^2 - 10^3 | Moderate; clinical possible with high target dwell time | Ibrutinib BTK ~5x10^3 (Pan 2007) |
| < 10^2 | Weak; needs warhead optimization | Reversibility likely dominates |
| <= 10 | Probably not covalent (or wrong residue) | Background rate vs GSH |

## Intrinsic Reactivity Assays

Before committing to a warhead, measure intrinsic reactivity (off-target risk):

```python
# Generic GSH stability assay readout - measure half-life of warhead with 10 mM GSH
# kinact_GSH from time-course of warhead disappearance
```

| Warhead | GSH t1/2 at 10 mM | Risk |
|---------|---------------------|------|
| Chloroacetamide | minutes | High (reacts with off-target Cys) |
| Acrylamide | hours | Moderate |
| Substituted acrylamide (alpha-Me) | days | Low (drug-like) |
| Nitrile | days | Low |
| Sulfonyl fluoride | hours-days | Variable |

The "**GSH-stable**" warhead (t1/2 > 4 hours) is the modern target for druglike TCIs.

## Covalent Docking Tools

| Tool | Approach | Strength | Fails when |
|------|----------|----------|------------|
| DOCKovalent (Wu 2014) | Constraint-based DOCK | Free, well-validated | Browser-based; small library |
| GOLD covalent (CCDC) | GOLD with covalent constraint | Commercial; selectivity | License cost |
| AutoDock 4 covalent | AD4 with covalent bond | Open source | Slower than Vina |
| CovDock (Schrödinger) | Glide-based + covalent | Commercial best | License cost |
| MOE covalent | Triposite Discovery | Commercial | License cost |
| HCovDock (Wang 2023) | Hierarchical fragment + covalent | Open; supports many residues | Newer, less validated |
| ICM-Pro covalent | Active site grid + covalent | Commercial; metal centers | License cost |

For open-source covalent docking, **HCovDock** (2023) is the modern alternative; **DOCKovalent** is the longstanding standard.

## Example: KRAS G12C Inhibitor Design Workflow

**Goal:** Decorate a co-crystal scaffold with a cysteine-targeting warhead and rank candidates by covalent efficiency.

**Approach:** Load scaffold SMILES, enumerate acrylamide-bearing analogs, filter by reactivity selectivity, dock under covalent constraint, and rank by kinact/Ki surrogates.

```python
from rdkit import Chem

# Step 1: scaffold from co-crystal (4LRW or AMG510)
scaffold_smi = 'c1ccc(C(=O)NC2=Nc3c(...)cnc23)cc1'
scaffold = Chem.MolFromSmiles(scaffold_smi)

# Step 2: enumerate analogs with acrylamide warhead
def add_acrylamide(scaffold, attachment_atom_idx):
    """Append acrylamide (-NC(=O)C=C) at a hydrogen position"""
    warhead = Chem.MolFromSmiles('NC(=O)C=C')
    # ... combine via Chem.RWMol or fragment combination
    pass

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
| AlphaFold3 / Boltz-2 binding pose | Geometric fit to reactive Cys |

**Goal:** Approximate Michael-acceptor reactivity from 2D structure without running DFT.

**Approach:** Parse the SMILES, locate the acrylamide substructure, count alpha-carbon substituents (more substitution lowers LUMO and slows reactivity), and return a negative count as a relative reactivity proxy.

```python
def acceptor_lumo_surrogate(smi):
    # Crude: count alpha-substituents to acrylamide; more substituents = lower reactivity
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    acryl_pat = Chem.MolFromSmarts('[CX3](=O)[CX3]=[CX3]')
    matches = mol.GetSubstructMatches(acryl_pat)
    if not matches:
        return None
    # Count substituents on alpha-C
    alpha_c = mol.GetAtomWithIdx(matches[0][1])
    n_subs = len([n for n in alpha_c.GetNeighbors() if n.GetIdx() not in matches[0]])
    return -n_subs  # crude proxy (lower = more reactive)
```

For real reactivity prediction, DFT calculations (LUMO energy, HOMO-LUMO gap) are needed.

## Per-Tool Failure Modes

### Wrong warhead for residue

**Trigger:** Acrylamide warhead targeted at Tyr.

**Mechanism:** Acrylamides are Cys-selective; do not form bonds with Tyr at physiological pH.

**Symptom:** No covalent adduct observed despite docking pose.

**Fix:** Match warhead to residue: acrylamide/chloroacetamide for Cys; sulfonyl fluoride for Lys/Tyr/Ser.

### Excessive reactivity (off-target)

**Trigger:** Chloroacetamide in drug-candidate context.

**Mechanism:** Too reactive; forms adducts with off-target Cys (e.g., GSH t1/2 < 30 min).

**Symptom:** Toxicity in cell-based assays; non-specific binding signal.

**Fix:** Replace with acrylamide (more selective); add alpha-substitution to acrylamide for tunable reactivity.

### Geometric mismatch

**Trigger:** Warhead positioned but Cys-Cβ distance > 6 Å.

**Mechanism:** Even with reactive warhead, geometric reach matters; Cys side chain has limited reach.

**Symptom:** No covalent labeling in mass spec despite predicted docking.

**Fix:** Validate by measuring distance from warhead to Cys-Cβ; redock with constrained covalent bond.

### Reversibility unintended

**Trigger:** Designed irreversible TCI but warhead is reversible.

**Mechanism:** Nitrile, aldehyde, boronate are reversible; equilibrium with non-covalent.

**Symptom:** Activity wanes after substrate washout in cellular assays.

**Fix:** Use truly irreversible warhead (acrylamide, chloroacetamide); or design for reversible covalent intentionally.

### kinact/Ki conflation

**Trigger:** Optimizing for IC50 instead of kinact/Ki.

**Mechanism:** Compounds with same IC50 differ in covalent efficiency.

**Symptom:** Apparently identical compounds have different in vivo PK.

**Fix:** Always measure kinact/Ki (kinetic assay); rank by covalent efficiency.

### DOCKovalent over-prediction

**Trigger:** Default DOCKovalent run.

**Mechanism:** Covalent constraint forces docking; many ligands "succeed" but are unrealistic.

**Symptom:** Many compounds pass docking; few label in vitro.

**Fix:** Post-filter by reactivity (chemoinformatics), geometric fit (Cys-Cβ distance), and PoseBusters.

## Reconciliation: Irreversible vs Reversible Covalent

| Aspect | Irreversible | Reversible covalent |
|--------|--------------|---------------------|
| Examples | KRAS G12C (acrylamide), BTK (ibrutinib) | Boronate (bortezomib), aldehyde (calpain inhibitors) |
| Toxicity profile | Off-target Cys labeling potential | Off-target equilibrium |
| Resistance mechanism | Mutation of reactive Cys | Mutation reduces affinity |
| Patent / IP | Stronger (specific bond) | Standard |
| When to choose | If Cys is hot-spot, conserved, druggable | If reversibility critical (e.g., proteasome) |

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Warhead not matching SMARTS | Different stereochemistry or charged | Use canonicalized + neutral SMARTS |
| DOCKovalent rejects ligand | No suitable Cys in pocket | Re-check residue accessibility |
| GSH adduct dominates | Warhead too reactive | Use less reactive warhead; or alpha-substitute |
| Off-target labeling in cells | Promiscuous warhead | Iterate warhead reactivity vs selectivity |
| Docking pose but no labeling | Geometric mismatch | Distance check; rotamer search |
| Reversible inhibitor not acting irreversibly | Wrong warhead class | Re-check reaction mechanism |
| HCovDock fails on PROTAC | Tool optimized for monomer covalent | Use specialized tools for bivalent |

## References

- Lonsdale & Ward, *Chem. Soc. Rev.* 47:3816 (2018) -- covalent revolution review.
- Singh et al., *Nat. Rev. Drug Discov.* 10:307 (2011) -- TCI design principles.
- Wu et al., *Bioinformatics* 30:i108 (2014) -- DOCKovalent.
- Wang et al., *Briefings Bioinform.* 24:bbac559 (2023) -- HCovDock.
- Cai et al., *J. Cheminformatics* 14:39 (2022) -- GOLD covalent toolkit.
- Backus et al., *Nat. Chem.* 8:530 (2016) -- proteome-wide covalent ABPP.
- Pettinger et al., *Angew. Chem. Int. Ed.* 56:15200 (2017) -- reactive warhead reactivity quantification.
- Schwartz et al., *Nat. Chem. Biol.* 10:1006 (2014) -- KRAS G12C disulfide-tethered fragments.

## Related Skills

- chemoinformatics/molecular-io - Parse warhead SMILES
- chemoinformatics/substructure-search - Warhead SMARTS detection
- chemoinformatics/virtual-screening - Pre-dock candidate non-covalent fit
- chemoinformatics/pose-validation - Validate covalent docking
- chemoinformatics/conformer-generation - Warhead conformer ensembles
- chemoinformatics/admet-prediction - ADMET of covalent leads
- chemoinformatics/molecular-descriptors - Reactivity surrogate descriptors

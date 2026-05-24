---
name: bio-protac-degraders
description: Designs PROTACs, molecular glues, and bivalent degraders with explicit handling of E3 ligase choice (VHL, CRBN, IAP, MDM2, KEAP1), linker design (length, composition, rigidity), ternary complex prediction (PRosettaC, DeepTernary, AlphaFold3 with constraints), cooperativity (alpha), DC50 / Dmax characterization, hook effect, and prediction-experiment reconciliation. Use when designing targeted protein degraders, planning linker SAR, predicting ternary complex stability, or building generative degrader workflows.
tool_type: python
primary_tool: PRosettaC
---

## Version Compatibility

Reference examples tested with: PRosettaC (web service), DeepTernary 1.0+, AlphaFold3 (constraints-enabled), Boltz-1 / Boltz-2, RDKit 2024.09+, OpenMM 8.1+ (for ternary MD).

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# PROTAC and Bivalent Degrader Design

Design bifunctional molecules (PROTACs) that recruit an E3 ubiquitin ligase to a target protein, inducing target ubiquitination and proteasomal degradation. PROTACs differ from traditional drugs: a stable **ternary complex** (target + PROTAC + E3) is required, not just target binding. The PROTAC field exploded post-2020 with clinical successes (ARV-471 estrogen-receptor degrader, ARV-110 androgen-receptor degrader). Postdoc-grade PROTAC design balances **target ligand binding**, **E3 ligand binding**, **linker geometry** (length, rigidity, chemistry), **cooperativity** (positive = ternary stable; negative = hook effect), and **cell permeability** (PROTACs are 800-1500 Da, often Lipinski-violating).

For target ligand design, see `chemoinformatics/virtual-screening` and `chemoinformatics/admet-prediction`. For linker-only enumeration, see `chemoinformatics/reaction-enumeration`. For generative linker design, see `chemoinformatics/generative-design`.

## E3 Ligase Choice

| E3 ligase | Ligand series | Best at | Limitations |
|-----------|---------------|---------|-------------|
| VHL | VL-269 (Salami 2018) | Surface-exposed targets | Tissue-restricted expression |
| CRBN (cereblon) | thalidomide, pomalidomide | Broad tissue expression | Off-target neosubstrates (IKZF1, SALL4) |
| IAP (XIAP, cIAP1) | SMAC mimetics (LCL161, etc.) | Apoptotic / IAP targets | Limited target scope |
| MDM2 | nutlin / idasanutlin | TP53 pathway | Limited target diversity |
| KEAP1 | DDB1-DCAF15-Keap1 | NRF2 pathway | Specialized use |
| RNF114 | EN450 | Newer; under exploration | Limited tooling |
| RNF4 | EN450-RNF4 | Newer; growing | Limited tooling |
| Cyclin E - SKP2 | New ligands emerging | Targeted cancers | Discovery-phase |

**Decision:** For first-generation PROTAC, **CRBN (cereblon)** is the most-developed (thalidomide-derived; broad tissue distribution). **VHL** is second-most-developed (more selective; tissue-restricted). For specialty targets, consider IAP / MDM2.

## Linker Design Principles

Linkers tune ternary complex geometry and stability:

| Property | Range | Effect |
|----------|-------|--------|
| Linker length | 8-30 atoms | Critical; geometry-dependent |
| Linker rigidity | Flexible (PEG) vs rigid (piperazine, pyridine) | Higher rigidity often reduces entropy penalty |
| Linker chemistry | PEG, alkyl, piperazine, triazole, ether, amide | PEG common; rigid for tighter binding |
| Click chemistry compatibility | Triazole compatible | Easy synthesis |
| MW range | PROTAC 800-1500 Da | Lipinski-violating but accepted |
| Polar atoms | 1-5 per linker | Permeability vs solubility balance |

**Critical:** The "Goldilocks linker length" is target-specific. Too short = ternary clash; too long = ternary entropy too high. Typically 12-20 atoms for surface-exposed targets.

## Decision Tree by Scenario

| Goal | E3 / linker | Tools |
|------|-------------|-------|
| First-generation PROTAC, surface-exposed target | CRBN + PEG linker (10-15 atoms) | PRosettaC for ternary prediction |
| Selective degrader (avoid off-target) | VHL + rigid linker | PRosettaC + cellular validation |
| BTK / IAP family targets | IAP-based PROTAC | Standard pipelines |
| Targeted protein degradation in cancer | CRBN or VHL | Standard clinical track |
| Novel target, no cryptic | Multiple E3 / linker variants | Combinatorial design + PRosettaC |
| Molecular glue (non-PROTAC) | CRBN-based | Distinct mechanism; pomalidomide-class |
| Optimize cooperativity (alpha > 1) | PRosettaC iterative; experimental | Ternary MD + ITC |
| Cell-active candidate | Standard development | PK + degradation cellular assays |

## Ternary Complex Prediction Tools

| Tool | Approach | Strength | Fails when |
|------|----------|----------|------------|
| PRosettaC | Rosetta + REMD on PROTAC | Mature; reliable | License Rosetta |
| DeepTernary | Equivariant deep learning | Fast; SE(3) | OOD chemistry |
| AlphaFold3 + constraints | Foundation model + restraints | High accuracy | Public access limited |
| Boltz-1 / Boltz-2 + constraints | Foundation model | Fast; flexible | Limited PROTAC training |
| HADDOCK | Restraint-based MD | Veteran | Manual restraint specification |
| Schrodinger Phase + Glide | Commercial | Production-ready | License cost |

**Decision:** For first-pass ternary modeling, **PRosettaC** is the standard. For prospective generative + ranking, **AlphaFold3** with constraints (if access) or **Boltz-1/2** (open). **DeepTernary** is a fast alternative when validated against benchmark.

## Cooperativity (Alpha)

Cooperativity quantifies how the ternary complex stabilizes (or destabilizes) the binary binding:

```
alpha = (Kd_binary,target) / (Kd_ternary,target)
```

- alpha > 1: positive cooperativity (ternary stronger than binary)
- alpha = 1: no cooperativity (independent binding)
- alpha < 1: negative cooperativity (mutual destabilization)

**Positive cooperativity** (alpha > 2-3) is the gold standard for PROTAC design. It means the ternary complex is stabilized by protein-protein interactions induced by the linker geometry.

Measure with ITC (isothermal titration calorimetry) or SPR/BLI titrations of binary vs ternary.

## DC50 / Dmax Characterization

In cellular assays:
- **DC50**: PROTAC concentration for 50% degradation (analogous to IC50)
- **Dmax**: maximum fraction degraded at any concentration

| Property | Good value | Notes |
|----------|------------|-------|
| DC50 | < 100 nM | Clinically meaningful |
| Dmax | > 80% | Sufficient depletion for phenotype |
| Hook effect | < 10x range | Linear dose-response over 100x |
| Hook concentration | typically > 10x DC50 | Above which binary complexes dominate |

**Hook effect**: at high PROTAC concentrations, binary complexes (PROTAC-target alone, PROTAC-E3 alone) dominate, and ternary complex formation drops. Dose-response curves are bell-shaped.

## Ternary Complex Modeling Workflow

**Goal:** Predict 3D structure of target-PROTAC-E3 ternary complex.

**Approach:**
1. Start with binary co-crystals: target + target-ligand pose; E3 + E3-ligand pose
2. Connect via linker enumeration (combinatorial)
3. Score by geometric feasibility (linker length, no clashes)
4. Refine with energy minimization

```python
# Pseudo-code workflow
def predict_ternary(target_pdb, target_ligand_sdf,
                    e3_pdb, e3_ligand_sdf, linker_smiles):
    # 1. Place binary complexes in same coordinate frame
    # 2. Enumerate linker connectivity from target-ligand exit vector to e3-ligand entry vector
    # 3. Score by total linker length, RMSD to expected geometry
    # 4. Minimize with restraint MD (5-20 ns)
    return ternary_poses
```

For production workflow, use PRosettaC web service or AlphaFold3 with chain-chain constraint inputs.

## Linker Length Calculation

```python
from rdkit import Chem
from rdkit.Chem import AllChem

def linker_distance(target_pose, e3_pose, target_ligand_smi, e3_ligand_smi, n_conf=20):
    """
    Compute distance between target-ligand exit vector and e3-ligand entry vector.
    Linker must span this distance with flexibility.
    """
    # Place target + e3 in same frame (from PRosettaC or HADDOCK)
    # Find target-ligand attachment atom and e3-ligand attachment atom
    # Compute distance
    pass

def required_linker_atoms(distance_A, rigidity='flexible'):
    """
    Estimate atoms needed in linker to span distance.
    Flexible (CC bond ~1.5 A) -> distance / 1.5
    Rigid (sp2 C-C ~1.4 A) -> distance / 1.4 + curvature factor
    """
    if rigidity == 'flexible':
        return int(distance_A / 1.4 * 1.2)  # 20% slack
    return int(distance_A / 1.4 * 1.5)
```

## Generative Linker Design

Combine REINVENT 4 linker mode with ternary prediction:

```toml
[run]
type = "linker"
prior_file = "priors/linker.prior"

[parameters]
input_smiles = "[*]NC(=O)..."  # Target ligand with attachment
linker_target_smiles = "[*]NC(=O)..."  # E3 ligand with attachment
n_steps = 200

[scoring_function.components]
type = "linker_length"
weight = 0.3
target_atoms = 14

type = "ternary_score"
weight = 0.5
ternary_predictor = "deepternary"  # API call

type = "qed"
weight = 0.2
```

REINVENT 4 supports linker generation; scoring needs ternary prediction (currently expensive). Modern best practice: combinatorial enumeration with PRosettaC scoring.

## Per-Tool Failure Modes

### PRosettaC -- inaccessible E3 in selected ligase

**Trigger:** Target's known binding mode incompatible with E3 ligase orientation.

**Mechanism:** PRosettaC samples linker conformations; can't bridge if E3 surface is occluded.

**Symptom:** Low ternary complex scores; high RMSD across replicates.

**Fix:** Try different E3 (CRBN vs VHL); use ARV-110 / ARV-471 family templates.

### DeepTernary -- novel chemotype

**Trigger:** Target ligand or E3 ligand outside training distribution.

**Mechanism:** DeepTernary learns from TernaryDB; novel chemotypes extrapolate.

**Symptom:** Predicted ternary complex unrealistic.

**Fix:** Validate against PRosettaC; use as fast screening, validate top-N with PRosettaC.

### Hook effect at low PROTAC concentration

**Trigger:** PROTAC linker too long; binary preferred energetically.

**Mechanism:** Negative cooperativity at low concentration; ternary unstable.

**Symptom:** Dose-response curve narrow; DC50 close to hook.

**Fix:** Shorter linker; increase cooperativity by adding protein-protein contact-promoting motifs.

### Insufficient cell permeability

**Trigger:** PROTAC > 1200 Da, polar, high TPSA.

**Mechanism:** Lipinski violations; permeability bottleneck.

**Symptom:** Cellular DC50 1000x worse than biochemical DC50.

**Fix:** Optimize linker for amphipathic profile; PROTACs often violate Ro5 but successful ones are < 1200 Da with < 130 TPSA.

### E3-target distance miscalculation

**Trigger:** Computing linker length from binary models without ternary refinement.

**Mechanism:** Binary structures don't capture ternary geometry.

**Symptom:** PROTACs synthesized at wrong linker length; no degradation.

**Fix:** Use PRosettaC-refined ternary distance; enumerate ±3 atoms around predicted optimal.

### Molecular glue vs PROTAC confusion

**Trigger:** Designing as PROTAC when target lacks defined ligand.

**Mechanism:** Molecular glues (like thalidomide for IKZF1) don't have target ligand component.

**Symptom:** Design too rigid; no degradation despite ternary prediction.

**Fix:** For targets without known ligand, consider molecular glue discovery instead.

## Reconciliation: PRosettaC vs AlphaFold3 (Constraint-Based)

| Aspect | PRosettaC | AlphaFold3 |
|--------|-----------|------------|
| Approach | REMD on linker + binary | Foundation model + restraints |
| Accuracy | High (validated benchmark) | Limited PROTAC data in training |
| Speed | Hours per ternary | Minutes |
| Open access | Yes (web service) | Public API, limited usage |
| Customization | Limited | Constraint flexibility |
| Decision | Standard for PROTAC | Promising; validate |

For PROTAC design, **PRosettaC + ITC/SPR experimental validation** is the current best practice.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| PRosettaC fails to converge | Linker too long or E3 surface buried | Try different E3; check binary structure quality |
| DeepTernary returns clashing pose | OOD chemotype | Use PRosettaC for refinement |
| AlphaFold3 ternary unrealistic | Constraint not specified properly | Use chain-chain distance constraints |
| PROTAC cell-active but biochemical degradation poor | Off-target neosubstrate (CRBN-IKZF) | Confirm target specificity via mass spec |
| Hook effect at low PROTAC concentration | Negative cooperativity | Optimize linker geometry |
| Synthesis too complex | Linker via 5+ step | Use click chemistry (triazole linker) |
| MW > 1500 Da issue | Excessive linker | Reduce linker length; rigid linker |

## References

- Békés et al., *Nat. Rev. Drug Discov.* 21:181 (2022) -- PROTAC clinical review.
- Drummond ML, Williams CI 2019 *J Chem Inf Model* 59(4):1634-1644 -- PROTAC ternary-complex modeling principles.
- Schapira et al., *Nat. Rev. Drug Discov.* 18:949 (2019) -- targeted protein degradation overview.
- Zaidman D, Prilusky J, Kessler N, Hadar D, Zarivach R, London N 2020 *J Chem Inf Model* 60(10):4894-4903 -- PRosettaC for ternary modeling (DOI 10.1021/acs.jcim.0c00589).
- Liu et al., *Nat. Commun.* (2025) -- DeepTernary deep learning.
- Saunders et al., *Sci. Rep.* (2025) -- PRosettaC vs AlphaFold3 comparison.
- Bondeson et al., *Nat. Chem. Biol.* 11:611 (2015) -- early PROTAC clinical readout.

## Related Skills

- chemoinformatics/molecular-io - Parse linker and ligand SMILES
- chemoinformatics/reaction-enumeration - Linker enumeration combinatorial
- chemoinformatics/generative-design - REINVENT linker mode
- chemoinformatics/conformer-generation - Ternary conformer sampling
- chemoinformatics/virtual-screening - Validate target ligand binding
- chemoinformatics/free-energy-calculations - Ternary ABFE / cooperativity
- chemoinformatics/admet-prediction - PROTAC ADMET specific challenges
- structural-biology/structure-io - PDB / mmCIF for ternary complex

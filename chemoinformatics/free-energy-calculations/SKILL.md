---
name: bio-free-energy-calculations
description: Performs alchemical free-energy calculations including relative binding free energy (RBFE / FEP+) and absolute binding free energy (ABFE) via OpenFE, FEP+, GROMACS, AMBER pmemd, and OpenMM with explicit lambda window scheduling, soft-core potentials, REST2 enhanced sampling, MBAR/BAR analysis, and cycle closure validation. Compares ML alternatives (Boltz-2 affinity, DeepDock). Use when ranking analogs by binding affinity beyond docking accuracy, performing prospective lead optimization, or validating SAR predictions.
tool_type: mixed
primary_tool: OpenFE
---

## Version Compatibility

Reference examples tested with: OpenFE 1.7+, OpenMM 8.1+, GROMACS 2024+, AMBER pmemd 22+, alchemlyb 2.1+, pymbar 4.0+, RDKit 2024.09+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `openfe --version`; `gmx --version`; `pmemd.cuda --version`

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Free Energy Calculations

Predict binding affinity differences (RBFE) or absolute binding affinities (ABFE) using alchemical free-energy methods. FEP+ (Schrödinger) is the commercial industry standard; OpenFE (Open Free Energy) is the open-source reference. Modern best practice achieves 1-2 kcal/mol RMSE vs experimental for well-set-up RBFE on rigid receptors. Boltz-2 affinity module (Wohlwend 2025) approaches FEP accuracy at 1000x speed on benchmarks, but FEP remains gold standard for production lead optimization.

For docking input poses, see `chemoinformatics/virtual-screening`. For pose validation before FEP, see `chemoinformatics/pose-validation`. For ML alternatives, see `chemoinformatics/ml-docking-rescoring`.

## FEP Method Taxonomy

| Method | Cost / pair | Accuracy | Use case | Fails when |
|--------|-------------|----------|----------|------------|
| FEP+ (Schrödinger) | hours-days GPU | 1-2 kcal/mol RMSE | Commercial lead opt | License cost |
| OpenFE RBFE | hours-days GPU | comparable to FEP+ | Open-source RBFE | Setup automation less mature |
| OpenFE ABFE | days GPU | 2-3 kcal/mol RMSE | Absolute affinity | Slower; more setup care |
| GROMACS RBFE | hours-days GPU | 1-2 kcal/mol | Power users, custom setup | Manual setup is error-prone |
| AMBER pmemd RBFE | hours-days GPU | 1-2 kcal/mol | Tradition; force-field maturity | Manual setup |
| FEP-SPell-ABFE | days GPU | 2-3 kcal/mol | Automated ABFE | Limited adoption |
| QligFEP v2.1 | minutes-hours | 1.5-3 kcal/mol | Q-based ligand FEP | Less standard |
| MM/PBSA | minutes | 3-5 kcal/mol RMSE | Endpoint, fast | Limited accuracy; entropy missing |
| MM/GBSA | minutes | 3-5 kcal/mol RMSE | Endpoint, faster than PBSA | Same caveats |
| Boltz-2 affinity | seconds GPU | 0.66 Pearson on FEP subset | ML alternative; 1000x faster | Novel chemotypes |
| ALEPB / EE-AMBER | days | 1-2 kcal/mol | Specialized | Limited tools |

**Decision:** For lead-optimization SAR validation, **OpenFE RBFE** (open) or **FEP+** (commercial) is the standard. For prospective discovery, MM/GBSA is a fast first-pass (3-5 kcal/mol RMSE); use FEP for top 10-50 candidates.

## Decision Tree by Scenario

| Scenario | Recommended workflow |
|----------|---------------------|
| Rank close analogs (R-group SAR) | RBFE via OpenFE (cycle: lig1↔lig2↔lig3) |
| Cross-scaffold ranking | ABFE per ligand; or coordinated RBFE with star network |
| Lead optimization 10-50 compounds | RBFE; perturbation-graph design |
| Single ligand affinity | ABFE (no reference needed) |
| Quick first-pass on top 1k | MM/GBSA after docking |
| Novel scaffold prospective | Boltz-2 affinity + FEP confirmation on top |
| Selectivity (target vs off-target) | RBFE on both proteins; report delta-delta-G |
| Allosteric vs orthosteric | ABFE comparable; check pose stability with MD |
| Ions / metal centers | Specialized force field (ZAFF, MCPB.py); not standard FEP |

## Relative Binding Free Energy (RBFE) Setup

**Goal:** Calculate delta-delta-G between two ligands (lig1 -> lig2) in pocket.

**Approach:** Alchemical transformation lig1 -> lig2 in both bound state (pocket + ligand + water) and unbound state (ligand + water alone). Thermodynamic cycle:

```
delta(delta-G_binding) = (delta-G_lig1->lig2 in pocket) - (delta-G_lig1->lig2 in solvent)
```

```python
# OpenFE simplified setup (real usage requires complete protocol setup)
from openfe import SmallMoleculeComponent, ProteinComponent, SolventComponent
from openfe.protocols.openmm_rfe import RelativeHybridTopologyProtocol

protein = ProteinComponent.from_pdb_file('receptor.pdb')
ligA = SmallMoleculeComponent.from_sdf_file('ligand_A.sdf')
ligB = SmallMoleculeComponent.from_sdf_file('ligand_B.sdf')
solvent = SolventComponent()

protocol = RelativeHybridTopologyProtocol(
    RelativeHybridTopologyProtocol.default_settings()
)
```

OpenFE's `setup` automates: mapping atoms between ligands (LOMAP), building hybrid topology, generating lambda windows, equilibration, production MD.

## Lambda Window Scheduling

| Stage | Lambda values | Purpose |
|-------|---------------|---------|
| Decoupling (vdW) | 0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0 | Turn off ligand vdW |
| Charging (Coulomb) | 0.0, 0.25, 0.5, 0.75, 1.0 | Turn off ligand partial charges |
| Restraint (ABFE only) | 0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0 | Boresch-style restraints |

Modern best practice uses 12-20 lambda windows per leg. Sampling at each window: 5-20 ns. Total simulation time per pair: 1-5 GPU-days.

## REST2 Enhanced Sampling

REST2 (Replica Exchange with Solute Tempering) is the de facto standard for FEP enhanced sampling. Scales solute-solute and solute-solvent interactions; allows ligand to overcome local minima.

In FEP+, REST2 region typically includes:
- The entire ligand
- Flexible binding-site loops
- Catalytic / ionic residues with high pKa shift potential

In OpenFE, REST2 is automatically applied to the alchemical region.

## MBAR/BAR Analysis

After production simulation, extract delta-G via MBAR (Multistate Bennett Acceptance Ratio) or BAR (Bennett Acceptance Ratio). MBAR uses data from all windows simultaneously; BAR uses adjacent windows.

```python
from alchemlyb.parsing import gmx
from alchemlyb.estimators import MBAR
import pandas as pd

u_nks = []
for window in range(12):
    df = gmx.extract_u_nk(f'window_{window}.xvg', T=300)
    u_nks.append(df)

u_nk = pd.concat(u_nks)
mbar = MBAR().fit(u_nk)
print(f'delta-G: {mbar.delta_f_.iloc[0, -1]:.2f} +/- {mbar.d_delta_f_.iloc[0, -1]:.2f} kcal/mol')
```

`alchemlyb` is the standard analysis package for FEP results from GROMACS, AMBER, OpenMM.

## Cycle Closure Analysis

Thermodynamic cycles must close (sum of edges = 0). Cycle closure error = root-mean-square error across closed cycles.

```python
def cycle_closure(rbfe_results, cycle):
    # cycle is list of edges, each (lig_i, lig_j, delta_g, sd)
    total = sum(d_g for _, _, d_g, _ in cycle)
    total_var = sum(sd**2 for _, _, _, sd in cycle)
    return total, total_var ** 0.5
```

Acceptable cycle closure: < 0.5 kcal/mol RMS. Higher indicates insufficient sampling or force-field issues.

## Absolute Binding Free Energy (ABFE)

ABFE computes delta-G of binding for a single ligand (no reference compound).

**Goal:** Compute Kd or Ki for a single ligand prospectively.

**Approach:** Decouple ligand from solvated state and from pocket-bound state separately; correction terms for analytical end states.

```bash
openfe absolute-free-energy run \
  --protein receptor.pdb \
  --ligand ligand.sdf \
  --output abfe_results/ \
  --n-lambda-charge 5 \
  --n-lambda-vdw 12 \
  --n-lambda-restraint 7
```

ABFE is harder than RBFE: requires Boresch-style restraints to keep ligand near pocket during decoupling. Restraint contribution must be analytically corrected.

ABFE cost: ~3x RBFE cost per ligand.

## MM/PBSA, MM/GBSA Endpoint Methods

Fast (<1 hour) alternative; lower accuracy (3-5 kcal/mol RMSE):

```bash
# MM/GBSA via AMBER MMPBSA.py
MMPBSA.py -i input.in -cp complex.parm7 -rp receptor.parm7 \
          -lp ligand.parm7 -y trajectory.nc
```

Sample input:
```
&general
  startframe = 100, endframe = 1000, interval = 10
&gb
  igb = 5
&pb
  istrng = 0.150
```

**Use case:** Rank order top 100 docking poses; MM/GBSA correlates ~0.5-0.7 with experimental binding; better than docking score (0.3-0.5) but worse than FEP (0.7-0.9).

## Force Field Selection

| Force field | Use for | Notes |
|-------------|---------|-------|
| OPLS4 (Schrödinger) | FEP+ default | Commercial; well-tested |
| OpenFF SAGE 2.x | OpenFE default | Open-source modern |
| GAFF2 | AMBER FEP | Use for ligand only; protein FF14SB |
| GAFF | Legacy | Replaced by GAFF2 |
| CGenFF | CHARMM-style FEP | CHARMM force-field family |
| ANI-2x | Mixed QM/MM | Experimental for FEP |
| MACE-OFF | Modern ML force field | Promising for FEP, limited tooling |

For OpenFE: defaults are SAGE 2.1.0 for ligand, FF14SB for protein, TIP3P for water. Override only if benchmarking.

## Per-Tool Failure Modes

### Insufficient sampling

**Trigger:** Lambda windows simulated < 5 ns each; ligand has slow rotamer change.

**Mechanism:** REST2 helps but isn't a panacea; some conformational changes take 100s of ns.

**Symptom:** Replicates disagree by > 1 kcal/mol; cycle closure > 1 kcal/mol RMS.

**Fix:** Increase per-window sampling to 10-20 ns; add REST2 to additional residues; check if pose is genuinely stable.

### Force-field artifacts

**Trigger:** Charged ligand or charged pocket residue.

**Mechanism:** GAFF2/SAGE may misparameterize unusual functional groups (perfluoro, charged sulfonate near Asp/Glu).

**Symptom:** Large RBFE error (>2 kcal/mol) for a specific transformation.

**Fix:** Visual inspection; check ligand topology with rdkit; consider non-bonded fix or fragment-specific parameters.

### Mapping ambiguity

**Trigger:** Two ligands differ in scaffold (not just R-groups).

**Mechanism:** LOMAP atom mapping may not find good correspondence; results from ambiguous mappings unreliable.

**Symptom:** Mapping score low; large dummy-atom count; cycle closure errors.

**Fix:** Manual mapping using OpenFE's editor; or use ABFE per ligand instead of RBFE.

### Restraint contribution wrong (ABFE)

**Trigger:** Boresch restraint applied to flexible region of ligand.

**Mechanism:** Analytical restraint correction assumes harmonic potential at well-defined minimum.

**Symptom:** ABFE differs by >3 kcal/mol from experiment systematically.

**Fix:** Choose Boresch restraint atoms from rigid ligand core; not flexible side chains.

### MM/GBSA -- bias from entropy missing

**Trigger:** Comparing ligands of very different size.

**Mechanism:** MM/GBSA misses entropy contribution; larger ligands appear more favorable.

**Symptom:** Larger ligands always rank higher.

**Fix:** Use MM/GBSA only for within-series ranking; supplement with FEP for cross-size.

### Boltz-2 affinity -- chemotype OOD

**Trigger:** Novel chemotype outside training distribution.

**Mechanism:** Boltz-2 is trained on PDBbind + ChEMBL; novel scaffolds extrapolate.

**Symptom:** Boltz-2 affinity and FEP affinity disagree.

**Fix:** Use Boltz-2 as a screen; validate top candidates with FEP. Treat Boltz-2 confidence band carefully.

## Reconciliation: FEP+ vs OpenFE

| Aspect | FEP+ | OpenFE |
|--------|------|--------|
| Force field | OPLS4 (proprietary) | SAGE 2.1.0 (open) |
| Workflow | Schrödinger GUI | Python CLI/API |
| Atom mapping | Automated LOMAP-style | LOMAP via openmm-tools |
| Reported accuracy | 1-2 kcal/mol RMSE | Comparable; emerging benchmarks |
| Cost | Schrödinger license | Free + compute time |
| Decision | Commercial team default | Open-source / academic / cost-sensitive |

For new groups, OpenFE 1.7+ is the recommended starting point; FEP+ is the gold standard for established pharma pipelines.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Lambda window simulation diverges | Bad initial pose | Re-relax pose with MM minimization first |
| Cycle closure > 1 kcal/mol | Insufficient sampling | Increase per-window time; check replicate convergence |
| MBAR returns NaN | Insufficient overlap between windows | Add intermediate lambda windows |
| Restraint contribution wrong | Boresch atoms on flexible region | Choose 3 atoms on rigid ligand core |
| GROMACS REST2 setup wrong | Hot region not specified | `-rest2-hot 'protein and resi 100-110'` style selection |
| ABFE under-estimates by ~3 kcal/mol | Forgetting analytical correction | Apply `delta-G_restraint_correction` term |
| MM/GBSA rmsd doesn't match docking | Different trajectory frames | Compute MM/GBSA on MD-relaxed pose |

## References

- Mey et al., *Living J. Comput. Mol. Sci.* 2:18378 (2020) -- alchemical free energy best practices.
- Wang et al., *J. Am. Chem. Soc.* 137:2695 (2015) -- FEP+ method.
- Open Free Energy (OpenFE) consortium 2023+ -- OpenFE framework. Cite the current release via the OpenFE Zenodo DOI (https://github.com/OpenFreeEnergy/openfe); the earlier "Henderson 2023 Comput Phys Commun" attribution could not be verified.
- Cournia et al., *J. Chem. Inf. Model.* 60:4153 (2020) -- RBFE for lead optimization.
- Aldeghi M et al -- ABFE protocols (consult current literature; the earlier "Aldeghi 2018 J Cheminform 10:43" citation could not be verified — Aldeghi's 2018 ABFE work appeared as a *Methods in Molecular Biology* book chapter).
- Wohlwend et al. (2025) -- Boltz-2 affinity prediction.
- Shirts & Chodera, *J. Chem. Phys.* 129:124105 (2008) -- MBAR.
- Bennett, *J. Comput. Phys.* 22:245 (1976) -- BAR.

## Related Skills

- chemoinformatics/virtual-screening - Source poses for FEP input
- chemoinformatics/pose-validation - PoseBusters-validate before FEP
- chemoinformatics/conformer-generation - Generate ligand 3D for FEP setup
- chemoinformatics/molecular-standardization - Standardize ligand before FEP
- chemoinformatics/ml-docking-rescoring - Boltz-2 affinity as alternative
- chemoinformatics/qsar-modeling - Surrogate models for high-throughput

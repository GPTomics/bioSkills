---
name: bio-free-energy-calculations
description: Performs alchemical free-energy calculations including relative binding free energy (RBFE / FEP+) and absolute binding free energy (ABFE) via OpenFE, FEP+, GROMACS, AMBER pmemd, and OpenMM with explicit lambda scheduling, soft-core potentials, MBAR/BAR analysis, cycle-closure validation, and protocol-appropriate enhanced sampling. Compares ML alternatives (Boltz-2 affinity, DeepDock). Use when ranking analogs by binding affinity beyond docking accuracy, performing prospective lead optimization, or validating SAR predictions.
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

Predict binding free-energy differences (RBFE) or standard binding free energies (ABFE) using alchemical methods. FEP+ is a commercial workflow and OpenFE is an open-source framework. Accuracy and cost vary substantially with system, perturbation, force field, setup, sampling, and evaluation design; report the protocol and benchmark relevant to the intended decision. The Boltz-2 report includes benchmark-specific comparisons with FEP methods but does not replace prospective validation on the project chemistry.

For docking input poses, see `chemoinformatics/virtual-screening`. For pose validation before FEP, see `chemoinformatics/pose-validation`. For ML alternatives, see `chemoinformatics/ml-docking-rescoring`.

## FEP Method Taxonomy

| Method | Cost / pair | Accuracy | Use case | Fails when |
|--------|-------------|----------|----------|------------|
| FEP+ (Schrödinger) | System- and protocol-dependent GPU cost | Published commercial RBFE workflow | Commercial lead optimization | License and reproducibility constraints |
| OpenFE RBFE | System- and protocol-dependent GPU cost | Open-source RBFE with documented protocols | Open-source campaigns | Mapping/setup/sampling require review |
| OpenFE ABFE | Generally more setup and sampling than one RBFE edge | Standard binding free energy | No congeneric reference ligand required | Restraints and end-state corrections |
| GROMACS / AMBER RBFE | Implementation-dependent | Custom alchemical workflows | Expert-controlled setup | Manual validation burden |
| FEP-SPell-ABFE | Protocol/system-dependent | Automated ABFE workflow | Evaluate published and project benchmarks | Limited adoption |
| QligFEP v2.1 | Protocol/system-dependent | Q-based ligand FEP | Evaluate published and project benchmarks | Different approximations/tooling |
| MM/PBSA / MM/GBSA | Lower-cost endpoint analysis | Approximate endpoint score | Exploratory within-series comparison | Entropy, sampling, and model dependence |
| Boltz-2 affinity | seconds GPU | 0.66 Pearson on reported FEP benchmark subset | ML alternative; reported >=1000x lower cost | Novel chemotypes |
| ALEPB / EE-AMBER | Protocol/system-dependent | Specialized methods | Evaluate matched evidence | Limited tooling |

**Decision:** For congeneric lead-optimization questions, evaluate a validated RBFE protocol and perturbation network. Use endpoint methods only for decisions supported by a project-specific benchmark. Candidate counts and escalation gates should follow compute budget, uncertainty, and prospective validation rather than a universal top-N rule.

## Decision Tree by Scenario

| Scenario | Recommended workflow |
|----------|---------------------|
| Rank close analogs (R-group SAR) | RBFE via OpenFE (cycle: lig1↔lig2↔lig3) |
| Cross-scaffold ranking | ABFE per ligand; or coordinated RBFE with star network |
| Congeneric lead-optimization set | RBFE with a connected, redundancy-aware perturbation graph |
| Single ligand affinity | ABFE (no reference needed) |
| Lower-cost exploratory ranking | A project-validated endpoint or ML method, followed by orthogonal confirmation |
| Novel scaffold prospective | Treat ML affinity as triage; validate selected decisions prospectively |
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

The protocol object does not itself choose an atom mapping or define a simulation. Create or inspect a mapping (Kartograf is the OpenFE 1.7 CLI default; LOMAP is also supported), construct bound and solvent `Transformation` objects, create their protocol DAGs, and run them through `openfe quickrun` or the documented Python execution interface. Always inspect the selected mapping before running.

## Lambda Window Scheduling

| Stage | Lambda values | Purpose |
|-------|---------------|---------|
| Decoupling (vdW) | 0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0 | Turn off ligand vdW |
| Charging (Coulomb) | 0.0, 0.25, 0.5, 0.75, 1.0 | Turn off ligand partial charges |
| Restraint (ABFE only) | 0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0 | Boresch-style restraints |

The 12-20 windows and 5-20 ns per-window ranges are repository starting ranges, not universal prescriptions. Select and extend them from overlap, exchange, and replicate-convergence diagnostics for the system; total cost therefore varies substantially.

## Enhanced Sampling

REST2 (Replica Exchange with Solute Tempering) is one enhanced-sampling approach used in some FEP workflows. It scales selected interactions to improve barrier crossing, but suitability and implementation are engine- and protocol-specific.

In FEP+, REST2 region typically includes:
- The entire ligand
- Flexible binding-site loops
- Catalytic / ionic residues with high pKa shift potential

FEP+ can use a configured REST2 region. OpenFE's `RelativeHybridTopologyProtocol` uses Hamiltonian replica exchange across its lambda states by default; that is not the same as REST2, and OpenFE does not automatically apply REST2. Use only enhanced-sampling modes supported and documented by the selected protocol and version.

## MBAR/BAR Analysis

After production simulation, extract delta-G via MBAR (Multistate Bennett Acceptance Ratio) or BAR (Bennett Acceptance Ratio). MBAR uses data from all windows simultaneously; BAR uses adjacent windows.

```python
from alchemlyb import concat
from alchemlyb.parsing import gmx
from alchemlyb.estimators import MBAR
from alchemlyb.postprocessors.units import to_kcalmol

u_nks = []
for window in range(12):
    df = gmx.extract_u_nk(f'window_{window}.xvg', T=300)
    u_nks.append(df)

u_nk = concat(u_nks)
mbar = MBAR().fit(u_nk)
delta_g = to_kcalmol(mbar.delta_f_).iloc[0, -1]
d_delta_g = to_kcalmol(mbar.d_delta_f_).iloc[0, -1]
print(f'delta-G: {delta_g:.2f} +/- {d_delta_g:.2f} kcal/mol')
```

`MBAR.delta_f_` and `d_delta_f_` are dimensionless (in kT) until explicitly converted. The parser shown above reads GROMACS XVG files. For other engines, use the engine-specific parser supported by the installed alchemlyb version.

## Cycle Closure Analysis

For a single directed thermodynamic cycle, the signed closure residual is the sum of its edges and should be consistent with zero within uncertainty. An RMS closure statistic requires residuals from multiple cycles and a stated aggregation convention.

```python
def cycle_closure_residual(cycle):
    # cycle is list of edges, each (lig_i, lig_j, delta_g, sd)
    total = sum(d_g for _, _, d_g, _ in cycle)
    total_var = sum(sd**2 for _, _, _, sd in cycle)
    return total, total_var ** 0.5
```

Interpret each closure residual relative to propagated edge uncertainties, replicate behavior, shared-edge correlations, network topology, and the decision supported. If reporting RMS across cycles, state which cycles were included and avoid treating correlated cycles as independent observations.

## Absolute Binding Free Energy (ABFE)

ABFE computes delta-G of binding for a single ligand (no reference compound).

**Goal:** Estimate the standard binding free energy of a single ligand prospectively. Conversion to an equilibrium dissociation constant requires an explicit standard-state convention; ABFE does not generically predict an assay `Ki`.

**Approach:** Decouple ligand from solvated state and from pocket-bound state separately; correction terms for analytical end states.

Use OpenFE's documented `AbsoluteBindingProtocol` workflow: construct the ligand and complex chemical systems, select and inspect the restraint setup, create the corresponding `Transformation` objects, serialize them with `Transformation.to_json()`, and execute each transformation with `openfe quickrun`. Do not substitute an ad hoc `absolute-free-energy` CLI; OpenFE does not provide that command.

ABFE is harder than RBFE: requires Boresch-style restraints to keep ligand near pocket during decoupling. Restraint contribution must be analytically corrected.

ABFE cost relative to RBFE depends on the protocols, number of legs/windows/repeats, and convergence requirements; estimate it from the explicit campaign plan.

## MM/PBSA, MM/GBSA Endpoint Methods

Lower-cost endpoint alternatives whose usefulness must be established on a matched benchmark:

```bash
# MM/GBSA via AMBER MMPBSA.py
MMPBSA.py -i input.in -cp complex.parm7 -rp receptor.parm7 \
          -lp ligand.parm7 -y trajectory.nc
```

Sample input:
```
&general
  startframe = 100, endframe = 1000, interval = 10
/
&gb
  igb = 5
/
&pb
  istrng = 0.150
/
```

**Use case:** Exploratory ranking when a matched retrospective benchmark shows the endpoint method supports the intended decision. Do not transfer generic correlation ranges across targets or protocols.

## Force Field Selection

| Force field | Use for | Notes |
|-------------|---------|-------|
| OPLS4 (Schrödinger) | FEP+ default | Commercial; well-tested |
| OpenFF 2.1.1 (Sage) | OpenFE 1.7 documented default | Inspect serialized settings; newer OpenFE releases use different defaults |
| GAFF2 | AMBER FEP | Use for ligand only; protein FF14SB |
| GAFF | Legacy | Replaced by GAFF2 |
| CGenFF | CHARMM-style FEP | CHARMM force-field family |
| ANI-2x | Mixed QM/MM | Experimental for FEP |
| MACE-OFF | Modern ML force field | Promising for FEP, limited tooling |

For OpenFE 1.7, the versioned documentation shows OpenFF 2.1.1 for the ligand and Amber-family protein/water XMLs including ff14SB and TIP3P. Inspect and serialize the actual protocol settings because defaults change between releases.

## Per-Tool Failure Modes

### Insufficient sampling

**Trigger:** Production length is insufficient for a slow ligand or protein degree of freedom.

**Mechanism:** Replica exchange improves state mixing but is not a panacea; some conformational changes remain slow.

**Symptom:** Replicates, time-sliced estimates, overlap/exchange diagnostics, or closure residuals are inconsistent with the reported uncertainty.

**Fix:** Increase sampling, inspect exchange and state overlap, run independent repeats, and investigate slow protein/ligand degrees of freedom. Use only protocol-supported enhanced sampling; check whether the pose is genuinely stable.

### Force-field artifacts

**Trigger:** Charged ligand or charged pocket residue.

**Mechanism:** GAFF2/SAGE may misparameterize unusual functional groups (perfluoro, charged sulfonate near Asp/Glu).

**Symptom:** A transformation is an outlier relative to experiment, replicates, or network consistency.

**Fix:** Visual inspection; check ligand topology with rdkit; consider non-bonded fix or fragment-specific parameters.

### Mapping ambiguity

**Trigger:** Two ligands differ in scaffold (not just R-groups).

**Mechanism:** LOMAP atom mapping may not find good correspondence; results from ambiguous mappings unreliable.

**Symptom:** Mapping score low; large dummy-atom count; cycle closure errors.

**Fix:** Manual mapping using OpenFE's editor; or use ABFE per ligand instead of RBFE.

### Restraint contribution wrong (ABFE)

**Trigger:** Boresch restraint applied to flexible region of ligand.

**Mechanism:** Analytical restraint correction assumes harmonic potential at well-defined minimum.

**Symptom:** ABFE shows a systematic offset or strong sensitivity to restraint choices.

**Fix:** Choose Boresch restraint atoms from rigid ligand core; not flexible side chains.

### MM/GBSA -- bias from entropy missing

**Trigger:** Comparing ligands of very different size.

**Mechanism:** MM/GBSA misses entropy contribution; larger ligands appear more favorable.

**Symptom:** Larger ligands always rank higher.

**Fix:** Use MM/GBSA only for within-series ranking; supplement with FEP for cross-size.

### Boltz-2 affinity -- chemotype OOD

**Trigger:** Novel chemotype outside training distribution.

**Mechanism:** Boltz-2 affinity training uses standardized public biochemical-assay data, including PubChem and ChEMBL sources, alongside its structural training. Novel targets, chemotypes, and assay contexts can still extrapolate.

**Symptom:** Boltz-2 affinity and FEP affinity disagree.

**Fix:** Use Boltz-2 as a benchmarked triage model and validate selected candidates with orthogonal computation and experiment. Do not interpret structure-confidence outputs as calibrated affinity intervals.

## Reconciliation: FEP+ vs OpenFE

| Aspect | FEP+ | OpenFE |
|--------|------|--------|
| Force field | OPLS4 (proprietary) | OpenFF 2.1.1 in versioned OpenFE 1.7 defaults; inspect serialized settings |
| Workflow | Schrödinger GUI | Python CLI/API |
| Atom mapping | Product workflow | Kartograf CLI default in OpenFE 1.7; LOMAP also supported |
| Reported accuracy | Benchmark-dependent | Benchmark-dependent; compare matched protocols and systems |
| Cost | Schrödinger license | Free + compute time |
| Decision | Commercial team default | Open-source / academic / cost-sensitive |

Choose OpenFE or a commercial workflow according to validated performance, auditability, available expertise, licensing, and integration requirements.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Lambda window simulation diverges | Bad initial pose | Re-relax pose with MM minimization first |
| Closure residual is inconsistent with propagated uncertainty or independent repeats | Sampling, mapping, force-field, or correlated-edge issue | Inspect signed residuals, overlap, mapping, and independent repeats before extending sampling |
| MBAR returns NaN | Insufficient overlap between windows | Add intermediate lambda windows |
| Restraint contribution wrong | Boresch atoms on flexible region | Choose 3 atoms on rigid ligand core |
| Slow binding-site rearrangement | Standard sampling does not cross the barrier | Increase sampling/repeats and use only engine- and protocol-documented enhanced sampling |
| ABFE systematic offset | Restraint, standard-state, sampling, or force-field issue | Inspect the protocol's documented restraint/free-energy terms and signs; do not invent an ad hoc correction variable |
| MM/GBSA rmsd doesn't match docking | Different trajectory frames | Compute MM/GBSA on MD-relaxed pose |

## References

- Mey ASJS et al., *Living J. Comput. Mol. Sci.* 2:18378 (2020) -- alchemical free-energy best practices (DOI 10.33011/livecoms.2.1.18378).
- Wang L et al., *J. Am. Chem. Soc.* 137:2695-2703 (2015) -- FEP+ method (DOI 10.1021/ja512751q).
- Open Free Energy developers. *OpenFE* software, Zenodo (2023-present) -- open-source alchemical free-energy framework (DOI 10.5281/zenodo.8344247).
- Cournia Z et al., *J. Chem. Inf. Model.* 60:4153-4169 (2020) -- rigorous ABFE as a final stage in virtual screening (DOI 10.1021/acs.jcim.0c00116).
- Aldeghi M, Bluck JP, Biggin PC. *Methods Mol. Biol.* 1762:199-232 (2018) -- beginner's guide to absolute alchemical ligand-binding free energies (DOI 10.1007/978-1-4939-7756-7_11).
- Passaro S et al. *bioRxiv* (2025) -- Boltz-2 affinity prediction preprint (DOI 10.1101/2025.06.14.659707).
- Shirts MR, Chodera JD. *J. Chem. Phys.* 129:124105 (2008) -- MBAR (DOI 10.1063/1.2978177).
- Bennett CH. *J. Comput. Phys.* 22:245-268 (1976) -- BAR (DOI 10.1016/0021-9991(76)90078-4).
- OpenFE 1.7 documentation: https://docs.openfree.energy/en/v1.7.0/
- alchemlyb documentation: https://alchemlyb.readthedocs.io/

## Related Skills

- chemoinformatics/virtual-screening - Source poses for FEP input
- chemoinformatics/pose-validation - PoseBusters-validate before FEP
- chemoinformatics/conformer-generation - Generate ligand 3D for FEP setup
- chemoinformatics/molecular-standardization - Standardize ligand before FEP
- chemoinformatics/ml-docking-rescoring - Boltz-2 affinity as alternative
- chemoinformatics/qsar-modeling - Surrogate models for high-throughput

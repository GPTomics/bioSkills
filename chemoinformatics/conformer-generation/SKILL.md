---
name: bio-conformer-generation
description: Generates 3D conformer ensembles using RDKit ETKDGv3 with knowledge-enhanced distance geometry, MMFF94/UFF force-field optimization, CREST + GFN2-xTB semi-empirical refinement, and macrocycle-aware torsion preferences. Provides explicit decision rules for single vs ensemble conformer use, RMSD pruning, energy windows, conformer count, and force-field choice. Use when preparing 3D ligands for docking, generating descriptor input for 3D QSAR, or sampling macrocycle/peptide conformational ensembles.
tool_type: mixed
primary_tool: RDKit
---

## Version Compatibility

Reference examples tested with: RDKit 2024.09+, xtb 6.7+, CREST 3.0+, OpenMM 8.1+ for follow-up MD.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `xtb --version`; `crest --version`

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Conformer Generation

Generate 3D conformer ensembles for molecules from 2D structures. The choice of method depends on molecule size, flexibility, and downstream use: ETKDG (Riniker & Landrum 2015) and its ETKDGv3 macrocycle update (Wang et al. 2020) are modern defaults for drug-like molecules, MMFF94/UFF provide fast energy minimization, and CREST + GFN2-xTB provide higher-cost semi-empirical sampling. A single conformer may be insufficient when the downstream result is conformation-sensitive; determine ensemble size by convergence of the downstream descriptor, alignment, or docking result.

For docking pose validation, see `chemoinformatics/pose-validation`. For free-energy methods (which require ensemble sampling), see `chemoinformatics/free-energy-calculations`.

## Conformer Method Taxonomy

| Method | Cost / mol | Quality | Use case | Fails when |
|--------|-----------|---------|----------|------------|
| ETKDGv3 + MMFF94 | Benchmark on actual molecules/hardware | Useful for many drug-like organics | Initial docking/descriptors | Difficult macrocycles, peptides, unsupported chemistry |
| ETKDGv3 + UFF | Fast | Different parameter coverage from MMFF94 | Fallback only after checking UFF parameters | Unsupported atom types; coordination chemistry |
| Omega (OpenEye) | Benchmark licensed workflow | Commercial conformer generator | Commercial pipelines | License cost and configured limits |
| Confab (Open Babel) | Benchmark on intended chemistry | Systematic torsion search | Alternative enumeration | Combinatorial growth and force-field dependence |
| RDKit ETKDGv3 + macrocycle preferences | Molecule-dependent | Macrocycle-aware embedding | Macrocyclic starting ensembles | Coverage remains molecule-dependent |
| CREST + GFN2-xTB | Molecule/settings-dependent | Semiempirical conformational sampling | Difficult flexible molecules | Computational cost; special chemistry |
| CREST + GFN-FF | Lower cost than GFN2-xTB | Force-field-level sampling | Exploratory sampling | Validate coverage and ordering for the chemistry |
| GeoMol (Ganea 2021) | Hardware/model-dependent | Learned conformer generation | Large-library research workflow | Training distribution and released-model coverage |
| TorsionNet (Gogineni 2020) | Hardware/model-dependent | Learned torsional search | Research workflow | Training distribution and implementation availability |
| MD sampling (OpenMM) | System/protocol-dependent | Dynamic sampling | Free energy, induced fit | Computational cost and convergence |

**Decision:** Start drug-like organic molecules with **ETKDGv3** and a parameter-checked MMFF94/MMFF94s optimization. Escalate difficult macrocycles, peptides, or highly flexible molecules to a validated CREST workflow when downstream convergence is inadequate. Benchmark ML generators on the intended chemistry before using them at scale.

## Decision Tree by Scenario

| Scenario | Starting method | Sampling and filtering decision |
|----------|-----------------|---------------------------------|
| Single initial 3D structure | ETKDGv3 + checked force field | Confirm embedding and minimization; downstream relaxation may still be required |
| Multi-conformer docking | ETKDGv3 ensemble | Increase sampling until pose recovery or enrichment is stable |
| 3D descriptors / pharmacophores | ETKDGv3 ensemble | Converge the reported statistic; justify energy/RMSD filters |
| Macrocycle / peptide | Macrocycle-aware ETKDG, then CREST if needed | Compare coverage against known conformers or downstream convergence |
| FEP input | Bound-pose-informed preparation and MD | Do not select solely by isolated-molecule conformer energy |
| Shape search | Query- and library-specific ensemble | Converge retrieval performance on a reference set |

## ETKDGv3 (Modern Default)

ETKDGv3 (Wang et al. 2020), building on the original ETKDG method (Riniker & Landrum 2015), incorporates experimental torsion preferences and updated macrocycle handling into distance geometry.

**Goal:** Generate an ensemble of 3D conformers from a SMILES with the modern default embedding algorithm.

**Approach:** Add explicit hydrogens, configure ETKDGv3 parameters (random seed, maximum embedding iterations, random coordinates), and embed multiple conformers via `EmbedMultipleConfs`.

```python
from rdkit import Chem
from rdkit.Chem import AllChem

def gen_conformers(smiles, n_conf=20, seed=42):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.useRandomCoords = True
    params.maxIterations = 1000
    ids = AllChem.EmbedMultipleConfs(mol, numConfs=n_conf, params=params)
    return mol, list(ids)
```

`useRandomCoords=True` can improve convergence for macrocycles and highly flexible molecules. On an `EmbedParameters` object the documented limit is `maxIterations`; `maxAttempts` is only present in legacy positional overloads.

## Force-Field Optimization

After embedding, minimize each conformer to a local minimum.

**Goal:** Reduce strain in each embedded conformer to a stable local minimum and record the resulting energies.

**Approach:** Build MMFF94s force-field parameters, minimize each conformer in place, and collect energies; fall back to UFF when MMFF94 cannot parameterize the molecule.

```python
def optimize_conformers(mol, conf_ids, force_field='mmff94'):
    results = []
    if force_field == 'mmff94':
        mmff_props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94s')
        if mmff_props is not None:
            for cid in conf_ids:
                ff = AllChem.MMFFGetMoleculeForceField(mol, mmff_props, confId=cid)
                if ff is None:
                    results.append({'conf_id': cid, 'status': 'force_field_failed'})
                    continue
                status = ff.Minimize(maxIts=1000)
                results.append({'conf_id': cid, 'energy': ff.CalcEnergy(),
                                'converged': status == 0, 'force_field': 'MMFF94s'})
            return results
        force_field = 'uff'
    if force_field == 'uff':
        if not AllChem.UFFHasAllMoleculeParams(mol):
            raise ValueError('Neither MMFF94s nor UFF covers this molecule')
        for cid in conf_ids:
            ff = AllChem.UFFGetMoleculeForceField(mol, confId=cid)
            if ff is None:
                results.append({'conf_id': cid, 'status': 'force_field_failed'})
                continue
            status = ff.Minimize(maxIts=1000)
            results.append({'conf_id': cid, 'energy': ff.CalcEnergy(),
                            'converged': status == 0, 'force_field': 'UFF'})
        return results
    raise ValueError(f'Unsupported force field: {force_field}')
```

**MMFF94 vs MMFF94s:** MMFF94s is the static-structure variant, with modified out-of-plane and torsional terms intended to preserve planarity in selected functional groups. It is useful for geometry optimization, but it is not a universally preferred replacement for MMFF94; record which variant was used and validate it for the downstream task.

**UFF (Universal Force Field):** UFF has different and often broader parameter coverage than MMFF94, but RDKit does not guarantee coverage for every molecule and generic UFF does not validate metal coordination chemistry. Call `UFFHasAllMoleculeParams` before use and use a chemistry-appropriate method for metal complexes.

## RMSD Pruning

Remove near-duplicate conformers within a chosen RMSD cutoff to keep the ensemble diverse:

```python
import numpy as np

def prune_conformers_rmsd(mol, conf_ids, rmsd_cutoff=0.5):
    n = len(conf_ids)
    keep = []
    for i, cid in enumerate(conf_ids):
        is_unique = True
        for kept_cid in keep:
            rmsd = AllChem.GetBestRMS(mol, mol, cid, kept_cid)
            if rmsd < rmsd_cutoff:
                is_unique = False
                break
        if is_unique:
            keep.append(cid)
    return keep
```

**Typical RMSD cutoff (Source / Rationale):**

| Cutoff | Use case | Source |
|--------|----------|--------|
| 0.5 Å | Drug-like ensemble for descriptors / docking | Repository clustering heuristic; validate for the downstream task |
| 1.0 Å | Drug-like ensemble for pharmacophore | Standard ROCS / pharmacophore practice |
| 1.5-2.0 Å | Macrocycles / peptides | Repository clustering heuristic for higher conformational freedom |
| 2.0+ Å | Cluster-centroid representative ensembles | Coarse representative sampling |

## Energy Window Filtering

Remove conformers above a project-justified energy cutoff only when the energy model and downstream purpose support that choice. Bound conformers can be strained relative to an isolated-molecule minimum.

```python
def filter_by_energy(mol, conf_ids, energies, window_kcal=10.0):
    min_e = min(energies)
    keep = []
    for cid, e in zip(conf_ids, energies):
        if e - min_e <= window_kcal:
            keep.append(cid)
    return keep
```

Treat any numerical energy window as a starting parameter. Calibrate it against conformer recovery or downstream metric convergence and record the energy method, solvent treatment, protonation state, and temperature assumptions.

## Macrocycle Handling

Macrocycles (>=12 atom rings) have distinct conformational issues: ETKDGv3 default knowledge base under-samples macrocycle torsions. Use macrocycle-specific torsion preferences:

```python
from rdkit.Chem import AllChem

def macrocycle_conformers(smiles, n_conf=200, seed=42):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f'Invalid SMILES: {smiles!r}')
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.useRandomCoords = True
    params.useMacrocycleTorsions = True
    params.useSmallRingTorsions = True
    params.maxIterations = 5000
    ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n_conf, params=params))
    if not ids:
        raise RuntimeError(f'No macrocycle conformers embedded for {smiles!r}')
    return mol, ids
```

For difficult macrocycles, CREST + GFN2-xTB is a useful higher-cost option; validate coverage against experimental or downstream evidence rather than treating one method as universally definitive.

## CREST + GFN2-xTB for High-Quality Sampling

CREST (Pracht et al. 2024) performs iterative meta-dynamics + GFN2-xTB optimization for conformer sampling.

**Goal:** Sample high-quality conformer ensembles for macrocycles, peptides, or molecules where ETKDGv3 + MMFF94 is inadequate.

**Approach:** Start from an RDKit-generated MMFF94-relaxed conformer, write to XYZ, and run CREST with GFN2-xTB driver to perform iterative meta-dynamics + reoptimization.

```bash
xtb mol.xyz --opt extreme
crest xtbopt.xyz --gfn2 --T 12 --ewin 6
```

**`--gfn2`**: use GFN2-xTB; validate its coverage and energy ordering for the chemistry of interest.
**`--gfn-ff`**: use GFN-FF; it can reduce computational cost, but benchmark its coverage and energy ordering for the intended molecules.
**`-ewin 6`**: 6 kcal/mol energy window above global min.
**`-T 12`**: use 12 CPU threads.

Output: `crest_conformers.xyz` with sampled ensemble.

**Workflow:** Start from RDKit ETKDGv3 + MMFF94 (cheap initial structure) -> save as XYZ -> CREST refinement.

```python
from rdkit import Chem
from rdkit.Chem import AllChem
import subprocess
from pathlib import Path

def crest_workflow(smiles, out_dir='crest_out'):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f'Invalid SMILES: {smiles!r}')
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.useRandomCoords = True
    if AllChem.EmbedMolecule(mol, params) == -1:
        raise RuntimeError(f'Initial 3D embedding failed for {smiles!r}')
    if not AllChem.MMFFHasAllMoleculeParams(mol):
        raise ValueError(f'MMFF parameters are unavailable for {smiles!r}')
    status = AllChem.MMFFOptimizeMolecule(
        mol, mmffVariant='MMFF94s', maxIters=1000
    )
    if status != 0:
        raise RuntimeError(f'Initial MMFF optimization did not converge (status {status})')

    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    input_path = out_dir / 'input.xyz'
    input_path.write_text(Chem.MolToXYZBlock(mol))
    # Because cwd is out_dir, pass the coordinate filename rather than out_dir/input.xyz.
    subprocess.run(['crest', input_path.name, '--gfn2', '-T', '12'],
                   cwd=out_dir, check=True)
    output_path = out_dir / 'crest_conformers.xyz'
    if not output_path.exists():
        raise FileNotFoundError(f'CREST did not produce {output_path}')
    return output_path
```

## Boltzmann Averaging of Properties

For ensemble descriptors (3D shape, dipole moment, polar surface area in 3D), Boltzmann-weight by energy:

```python
import numpy as np

def boltzmann_weights(energies, T=300.0):
    energies = np.array(energies)
    kt = 0.001987 * T  # kcal/mol at 300K
    rel = energies - energies.min()
    w = np.exp(-rel / kt)
    return w / w.sum()

def boltzmann_average(values, energies, T=300.0):
    w = boltzmann_weights(energies, T)
    return float(np.sum(np.array(values) * w))
```

Boltzmann populations require energies that approximate the relevant thermodynamic state, including conformer degeneracy and environmental effects where important. Raw gas-phase MMFF/UFF minima are exploratory surrogates, not generally validated populations.

## ML-Based Conformer Generation and Search

Research methods such as GeoMol generate molecular conformations directly, while TorsionNet uses reinforcement learning to search torsional conformational space. They are distinct approaches and should not be presented as interchangeable drop-in generators.

```python
# Pseudo-code for GeoMol-style ML conformer generation
# (Requires pre-trained model + dependencies)
# from geomol import generate_conformers
# conformers = generate_conformers(smiles, n_conformers=10)
```

**Trade-off:** Performance and coverage depend on the released model, training distribution, conformer definition, and benchmark. Verify the maintained implementation and benchmark it against ETKDGv3 or CREST on the intended chemistry before using it operationally; do not infer broad macrocycle or organometallic coverage from drug-like benchmarks.

## Per-Tool Failure Modes

### ETKDGv3 -- failed embedding

**Trigger:** Macrocycle, highly constrained polycyclic, or sterically crowded molecule.

**Mechanism:** Distance geometry cannot find a consistent 3D structure within the configured embedding iterations.

**Symptom:** `EmbedMolecule` returns -1; `EmbedMultipleConfs` returns empty list.

**Fix:** Set `useRandomCoords=True`, increase `maxIterations`; for macrocycles, set `useMacrocycleTorsions=True`. As fallback, use CREST.

### MMFF94 -- parameter missing

**Trigger:** Molecule contains element not parameterized (transition metals, certain S+ species).

**Mechanism:** MMFF94 only covers H, C, N, O, F, Si, P, S, Cl, Br, I + select cations.

**Symptom:** `MMFFGetMoleculeProperties` returns None; optimization silently no-ops.

**Fix:** Fall back to UFF; or for metals, use GFN2-xTB.

### Conformer ensemble too small

**Trigger:** `n_conf=10` for a flexible molecule (>5 rotatable bonds).

**Mechanism:** A fixed small ensemble can miss relevant minima for flexible molecules.

**Symptom:** New conformers continue to change cluster populations or the downstream result.

**Fix:** As a repository starting heuristic, use n_conf = max(10, 5 * NumRotatableBonds + 10), then increase sampling until the ensemble is stable for the downstream metric.

### Single-conformer 3D descriptor

**Trigger:** Calculating 3D descriptors from a single conformer.

**Mechanism:** Some 3D descriptors vary materially across conformers.

**Symptom:** Same molecule produces different 3D descriptors on rerun.

**Fix:** Converge the descriptor over an ensemble and report the chosen summary. Use Boltzmann weighting only with a justified population model.

### CREST -- timeout on flexible molecule

**Trigger:** Cyclosporin or large peptide.

**Mechanism:** CREST metadynamics scales poorly with rotational complexity.

**Symptom:** Hours of CPU time per molecule; incomplete sampling.

**Fix:** Use `--gfn-ff` for an exploratory lower-cost run or shorten metadynamics with the documented `--mdlen`/`--len` option. `--noopt` disables input pre-optimization; it does not skip metadynamics.

### GFN2-xTB conformer reordering

**Trigger:** Comparing conformer energies between GFN2-xTB and DFT.

**Mechanism:** GFN2-xTB is parameterized for energies; relative conformer ordering can differ from DFT by 1-2 kcal/mol.

**Symptom:** "Wrong" conformer reported as global minimum vs DFT reference.

**Fix:** For high-stakes work, re-rank top GFN2-xTB conformers with DFT single-points (e.g., r2SCAN-3c).

## Reconciliation: ETKDGv3 vs CREST

| Use case | ETKDGv3 | CREST |
|----------|---------|-------|
| Drug-like organic molecule | Efficient starting point | Higher-cost comparison when convergence fails |
| Highly flexible molecule | Increase and convergence-test sampling | Useful alternative sampling strategy |
| Macrocycle or peptide | Try macrocycle-aware settings and validate | Often useful, but not automatically sufficient |
| Population-weighted descriptors | Requires justified energy/population model | Higher-level energy still requires thermodynamic validation |
| FEP input | Useful for initial coordinates | Does not replace bound-pose and MD preparation |

For ETKDGv3 ensembles, compare a representative subset with an orthogonal method or experimental conformers and judge adequacy using the downstream metric; no single cross-method RMSD establishes completeness.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| `EmbedMolecule` returns -1 | Embed failed | Set `useRandomCoords=True`; raise `params.maxIterations` |
| MMFFOptimize no-op | MMFF parameters missing | Use UFF fallback |
| All conformers identical | Stiff molecule | OK; molecule is rigid |
| Conformers physically wrong | Stereochemistry lost | Re-add explicit stereo before embedding |
| 3D descriptors differ per run | Random seed not set | `params.randomSeed = 42` |
| CREST out-of-memory | Search/ensemble too large | Reduce `--T`, shorten documented sampling, or lower `--ewin` to retain fewer structures |
| Macrocycle ring inverted | Default torsion preferences wrong | Set `useMacrocycleTorsions=True` |
| AddHs not called | Implicit H not embedded | `mol = Chem.AddHs(mol)` before EmbedMolecule |

## References

- Hawkins et al., *J. Chem. Inf. Model.* 50:572-584 (2010) -- OMEGA conformer sampling (DOI 10.1021/ci100031x).
- Riniker & Landrum, *J. Chem. Inf. Model.* 55:2562-2574 (2015) -- original ETKDG method (DOI 10.1021/acs.jcim.5b00654).
- Wang S, Witek J, Landrum GA, Riniker S. *J. Chem. Inf. Model.* 60:2044-2058 (2020) -- ETKDGv3 macrocycle update (DOI 10.1021/acs.jcim.0c00025).
- Halgren TA, *J. Comput. Chem.* 17:490-519 (1996) -- MMFF94 force field (DOI 10.1002/(SICI)1096-987X(199604)17:5/6%3C490::AID-JCC1%3E3.0.CO;2-P).
- Rappe AK et al., *J. Am. Chem. Soc.* 114:10024-10035 (1992) -- UFF (DOI 10.1021/ja00051a040).
- Pracht P et al., *J. Chem. Phys.* 160:114110 (2024) -- CREST 3.0 (DOI 10.1063/5.0197592).
- Bannwarth C, Ehlert S, Grimme S. *J. Chem. Theory Comput.* 15:1652-1671 (2019) -- GFN2-xTB (DOI 10.1021/acs.jctc.8b01176).
- RDKit force-field API: https://www.rdkit.org/docs/source/rdkit.Chem.rdForceFieldHelpers.html
- CREST command-line documentation: https://crest-lab.github.io/crest-docs/page/documentation/keywords.html
- xTB documentation: https://xtb-docs.readthedocs.io/
- Ganea et al., *NeurIPS* (2021) -- GeoMol ML conformer generation.
- Gogineni T et al., *NeurIPS* 33 (2020) -- TorsionNet learned conformational search.

## Related Skills

- chemoinformatics/molecular-io - Parse molecules
- chemoinformatics/molecular-standardization - Standardize before embedding
- chemoinformatics/molecular-descriptors - 3D descriptors from ensembles
- chemoinformatics/shape-similarity - Multi-conformer 3D shape matching
- chemoinformatics/virtual-screening - Generate 3D ligands for docking
- chemoinformatics/free-energy-calculations - Sample conformers for MD setup
- chemoinformatics/pharmacophore-modeling - 3D pharmacophore from ensembles

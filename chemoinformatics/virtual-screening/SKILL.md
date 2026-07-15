---
name: bio-virtual-screening
description: Performs structure-based virtual screening using AutoDock Vina, SMINA, GNINA (CNN scoring), and DiffDock-L hybrid workflows with explicit choice rules across rigid vs flexible docking, cross-docking vs self-docking, binding-site detection (P2Rank, fpocket), receptor preparation (PDB2PQR, PROPKA), ligand preparation (meeko, OpenBabel), and ultralarge-library screening (ZINC22, Enamine REAL). Use when screening chemical libraries against a protein target to find candidate binders, ranking docking poses, or selecting a docking workflow for a specific scenario.
tool_type: python
primary_tool: AutoDock Vina
---

## Version Compatibility

Reference examples tested with: AutoDock Vina 1.2.5+, SMINA 2020-12+, GNINA 1.1+ for `rescore` (GNINA 1.3+ for the six-mode interface documented below), RDKit 2024.09+, meeko 0.5+, P2Rank 2.4+, ProDy 2.4+, pdb2pqr 3.6+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `vina --version`; `gnina --version`; `smina --version`

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Virtual Screening

Screen chemical libraries against protein targets via molecular docking. Vina is the de-facto default, SMINA adds flexibility (Vinardo scoring, custom scoring), and GNINA adds CNN-based pose scoring (Top-1 redock 58%->73% over Vina, cross-dock 27%->37%). Deep-learning docking (DiffDock-L, EquiBind, NeuralPLexer) competes in pose accuracy, but physical validity is method- and dataset-dependent; the workflow therefore combines ML pose sampling with classical scoring and explicit geometry checks. For ultralarge libraries (>1M), library preparation, hierarchical filtering, and HPC orchestration become the limiting steps.

For pose physical-validity QC, see `chemoinformatics/pose-validation`. For ML-driven docking + rescoring, see `chemoinformatics/ml-docking-rescoring`. For covalent docking, see `chemoinformatics/covalent-design`. For affinity calculations (FEP), see `chemoinformatics/free-energy-calculations`.

## Docking Tool Taxonomy

| Tool | Scoring | Speed (sec/lig) | Best at | Fails when |
|------|---------|-----------------|---------|------------|
| AutoDock Vina 1.2 | Vina (empirical) | Hardware- and settings-dependent | Open, well-characterized baseline | Cross-dock; cryptic pockets; metal centers |
| SMINA | Vina + flexible + custom | Hardware- and settings-dependent | Custom scoring; flexible side chains | Same Vina-scoring caveats |
| Vinardo | Modified Vina scoring | Hardware- and settings-dependent | Alternative empirical score | Validate on target-relevant controls |
| GNINA 1.1 | CNN or Vina scoring | GPU- and settings-dependent | CNN-assisted pose ranking | Validate transfer to the target and chemotype |
| AutoDock 4 | AD4 + grid maps | Hardware- and settings-dependent | Legacy reference | More setup than Vina |
| DOCK 6/7 | DOCK + Amber | Hardware- and settings-dependent | UCSF DOCK ecosystem | Steep learning curve |
| Glide (Schrodinger) | GlideScore | License and hardware-dependent | Commercial docking workflow | License cost |
| GOLD (CCDC) | GOLDScore / ChemScore | License and hardware-dependent | Commercial workflow; metal options | License cost |
| FlexX (BioSolveIT) | FlexX | License and hardware-dependent | Fragment-based placement | License cost |
| rDock | rDock | Hardware- and settings-dependent | Open-source alternative | Validate maintenance and target fit |
| DiffDock-L | Diffusion-generative | GPU- and settings-dependent | Pose sampling for cross-docking | Validate geometry with PoseBusters; see ml-docking-rescoring |
| EquiBind | Equivariant NN | GPU- and settings-dependent | Single-shot pose generation | Requires independent geometry and ranking checks |
| Boltz-2 + GNINA rescore | Foundation model + CNN | GPU- and settings-dependent | Experimental multi-model workflow | Benchmark each evidence stream independently |

**Decision:** Use Vina as an open baseline and consider GNINA CNN rescoring when target-relevant redocking or cross-docking controls support it. For large libraries, calibrate a hierarchical Vina -> GNINA -> higher-cost follow-up workflow on measured enrichment, throughput, and retained chemotype diversity.

## Decision Tree by Scenario

| Scenario | Recommended workflow |
|----------|---------------------|
| Self-dock against known ligand pocket | GNINA `gnina --cnn_scoring rescore` |
| Cross-dock to apo or related-target structure | DiffDock-L pose + GNINA rescore + PoseBusters |
| Ultralarge library (10M+) | Calibrated hierarchical screen: property/alert triage -> Vina -> measured top fraction to GNINA -> higher-cost follow-up |
| Cryptic pocket / induced fit | Receptor-ensemble docking and, where appropriate, a separately validated complex-prediction model |
| Allosteric / undefined site | P2Rank for pocket detection -> ensemble dock all pockets |
| Metal-coordinated ligand | GOLD (commercial) or manually parameterize Vina metal scoring |
| Covalent inhibitor | See `chemoinformatics/covalent-design`: DOCKovalent, HCovDock |
| Fragment screen (<300 Da) | rDock or constrained Vina with seed atoms |
| Hit-to-lead refinement | Use co-crystal structure if available; MD-relaxed receptor; FEP for affinity |

## Receptor Preparation

**Goal:** Convert a protein PDB into a docking-ready format with correct protonation, missing atoms, and removed waters.

**Approach:** Decide which ligands, cofactors, metals, and structural waters to retain -> fill missing heavy atoms with a structure-repair tool such as PDBFixer -> use PROPKA/PDB2PQR plus manual review to assign pH-dependent protonation -> assign the charge model required by the docking workflow -> prepare receptor PDBQT with a documented AutoDock-compatible tool.

```python
import subprocess
from pathlib import Path

def prepare_receptor(repaired_pdb, pdbqt_out, pH=7.4):
    # Decide which waters/cofactors/metals to retain before this function.
    base = str(Path(repaired_pdb).with_suffix(''))
    pqr_file = f'{base}_pH{pH}.pqr'
    subprocess.run(['pdb2pqr', '--ff=AMBER', f'--with-ph={pH}',
                    repaired_pdb, pqr_file], check=True)
    output_basename = str(Path(pdbqt_out).with_suffix(''))
    subprocess.run(['mk_prepare_receptor.py', '--read_pqr', pqr_file,
                    '-o', output_basename, '-p'], check=True)
    return pdbqt_out
```

**Common pitfall:** Forgetting to add hydrogens at protein pH (7.4) but using pH 7.0 ligand charges. Hist mistakenly protonated. Use PROPKA + manual review of catalytic residues.

## Ligand Preparation

**Goal:** Generate a 3D, docking-ready ligand file from SMILES with appropriate protonation and conformation.

**Approach:** Supply a documented protomer/tautomer state generated by an appropriate pKa/protomer workflow -> parse it with RDKit -> embed 3D with ETKDGv3 -> minimize with MMFF94 -> write PDBQT with Meeko. `MolFromSmiles` parses the supplied state and `Uncharger` neutralizes formal charges; neither predicts protonation at pH 7.4.

```python
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

def prepare_ligand(smiles, pdbqt_out):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f'invalid SMILES: {smiles}')
    mol = Chem.AddHs(mol)
    embed_status = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    if embed_status != 0:
        raise RuntimeError('ETKDGv3 failed to generate a ligand conformer')
    if not AllChem.MMFFHasAllMoleculeParams(mol):
        raise ValueError('MMFF94 parameters are unavailable for this ligand')
    optimization_status = AllChem.MMFFOptimizeMolecule(mol)
    if optimization_status != 0:
        raise RuntimeError('MMFF94 ligand optimization did not converge')

    # meeko 0.5+ API: prepare() returns a list of MoleculeSetup objects;
    # use PDBQTWriterLegacy.write_string() to materialize the PDBQT block.
    mk_prep = MoleculePreparation()
    setups = mk_prep.prepare(mol)
    pdbqt_text, is_ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not is_ok:
        raise RuntimeError(f'meeko PDBQT export failed: {err}')
    with open(pdbqt_out, 'w') as f:
        f.write(pdbqt_text)
    return pdbqt_out
```

`meeko` (AutoDock developers' tool) handles torsion tree creation, rotamer flagging, and PDBQT writing -- preferred over Open Babel's PDBQT writer. Note: meeko 0.5+ separated the writer (`PDBQTWriterLegacy`) from `MoleculePreparation`; older code using `prep.write_pdbqt_file()` is deprecated.

## Binding Site Detection

When the binding pocket is not known (apo target, novel allosteric site):

| Tool | Approach | Output |
|------|----------|--------|
| P2Rank (Krivak 2018) | ML on protein surface descriptors | Ranked pocket list with center coords |
| fpocket (Le Guilloux 2009) | Voronoi tessellation | Pocket descriptor list |
| DoGSiteScorer | Geometric + drugability | Pocket list with score |
| AutoSite (Vina) | Affinity map clustering | Pocket centers |
| AlphaFill | Transplant ligands/cofactors from homologous experimental structures into AlphaFold models | Plausible binding-site components for review |

```bash
prank predict -f receptor.pdb -o pockets/
```

P2Rank output `<receptor>_predictions.csv` lists pocket centers with scores. The highest model score does not identify a pocket as orthosteric or biologically relevant; verify ranked pockets against co-crystal, mutagenesis, SAR, or other structural evidence.

## Vina Docking (Single Ligand)

```python
# AutoDock Vina Python API requires Vina 1.2+; for Vina 1.1 use subprocess CLI:
# subprocess.run(['vina', '--receptor', ..., '--ligand', ..., '--center_x', ...], check=True)
from vina import Vina

def dock_single(receptor_pdbqt, ligand_pdbqt, center, box_size,
                exhaustiveness=8, n_poses=10):
    v = Vina(sf_name='vina')
    v.set_receptor(receptor_pdbqt)
    v.set_ligand_from_file(ligand_pdbqt)
    v.compute_vina_maps(center=center, box_size=box_size)
    v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
    return v.energies(), v.poses()
```

**Exhaustiveness:** `8` is the Vina default. Increasing it increases search effort, but runtime and pose recovery depend on hardware, ligand flexibility, box size, and software version. Benchmark settings such as 8, 16, 32, and 64 on target-relevant controls instead of assigning universal timing or quality labels.

Vina's `rmsd_lb` and `rmsd_ub` are lower and upper heavy-atom RMSD bounds between a reported mode and the best-scoring mode; the bounds differ in how symmetry-equivalent atoms are handled. They are not pose-versus-experimental-reference RMSDs. Use an external symmetry-aware RMSD to a reference pose for accuracy QC.

## GNINA with CNN Scoring (modern default)

```bash
gnina -r receptor.pdb -l ligand.sdf \
      --autobox_ligand reference_ligand.sdf \
      --cnn_scoring rescore \
      -o poses.sdf.gz \
      --num_modes 9 --exhaustiveness 8
```

`--cnn_scoring`:
- `none`: no CNN; use the selected empirical scoring function throughout
- `rescore` (default): use empirical scoring during the search, then CNN-rerank the final poses; least computationally expensive CNN option
- `refinement`: use the CNN to refine poses after Monte Carlo chains and to rank the final poses; approximately 10 times slower than `rescore` on a GPU in the official documentation
- `metrorescore`: use CNN scoring in the Metropolis search and rescore the resulting poses
- `metrorefine`: use CNN scoring in the Metropolis search and refine the resulting poses
- `all`: use the CNN scoring function throughout; the official documentation describes this as extremely computationally intensive and not recommended

The six choices above are from GNINA 1.3. Earlier releases expose a smaller set; check `gnina --help` for the installed executable rather than assuming every mode is available.

`--autobox_ligand`: define box from reference ligand SDF/PDB. Otherwise specify `--center_x/y/z` + `--size_x/y/z`.

**Critical:** GNINA distributions include multiple named CNN models/ensembles rather than one universally described "PDBbind 2019" model. Record the selected model or ensemble and validate it with known co-crystal redocking and, when relevant, cross-docking controls.

## Virtual Screening Pipeline (Hierarchical)

**Goal:** Screen 10M-compound library down to top-1k candidates for follow-up.

**Approach:** Three-stage filter. The 1% and top-1000 selections below are repository starting heuristics; choose production cutoffs from target-relevant enrichment, diversity, and throughput measurements.

Pseudo-code skeleton (orchestrator). Each helper function delegates to a dedicated skill: drug-likeness filter to `chemoinformatics/admet-prediction`, single-ligand Vina/GNINA to `dock_single` defined earlier in this skill, PoseBusters QC to `chemoinformatics/pose-validation`.

```python
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
from functools import partial

# Stub helpers to be implemented per project; see the cross-referenced skills.
def drug_like_filter(df):
    raise NotImplementedError('Implement via chemoinformatics/admet-prediction (Lipinski+Veber+PAINS)')
def vina_dock(smi, receptor_pdbqt, center, box):
    raise NotImplementedError('Wrap dock_single() above; return best affinity')
def gnina_rescore(smi, receptor_pdbqt, center, box):
    raise NotImplementedError('Wrap gnina --cnn_scoring rescore subprocess call')
def pose_validate(df):
    raise NotImplementedError('Implement via chemoinformatics/pose-validation (PoseBusters)')

def vs_pipeline(library_smi, receptor_pdbqt, center, box, output_dir, n_workers=16):
    df = pd.read_csv(library_smi)
    df_stage1 = drug_like_filter(df)

    worker = partial(vina_dock, receptor_pdbqt=receptor_pdbqt,
                     center=center, box=box)
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        affinities = list(ex.map(worker, df_stage1['smiles']))
    df_stage1['vina_affinity'] = affinities
    df_stage2 = df_stage1.nsmallest(int(len(df_stage1) * 0.01), 'vina_affinity')

    df_stage2['gnina_affinity'] = df_stage2['smiles'].apply(
        lambda smi: gnina_rescore(smi, receptor_pdbqt, center, box))
    df_stage3 = df_stage2.nsmallest(1000, 'gnina_affinity')

    return pose_validate(df_stage3)
```

For very large libraries, use a restartable scheduler-backed workflow and measure throughput on a representative tranche. Record hardware, software version, box dimensions, ligand flexibility, and failure rate with every throughput estimate.

## Ultralarge Library Screening (ZINC22, Enamine REAL)

| Library | Scope | Typical access | Verification requirement |
|---------|-------|----------------|--------------------------|
| ZINC22 | Purchasable and make-on-demand compounds | Tranche/download interfaces | Record the tranche query and retrieval date |
| Enamine REAL | Make-on-demand compounds | Provider files or search interface | Record product-space release and retrieval date |
| Enamine HTS | Screening collection | Provider files | Confirm current stock/version with the provider |
| Mcule | Aggregated purchasable compounds | Provider search/export | Record filters and retrieval date |
| ChEMBL | Curated compounds and bioactivities | Versioned database release | Record ChEMBL release and extraction query |

Library sizes and availability change frequently. Obtain counts from the provider or versioned database at execution time rather than copying a static total into a workflow.

For ultralarge VS, the following percentages and thresholds are repository starting heuristics that must be calibrated for the target and library:
1. Apply a documented property/alert policy while retaining flagged and rejected counts
2. If known actives exist, test a permissive 2D-similarity prefilter such as ECFP4 Tanimoto >=0.4 and measure active/chemotype retention
3. Vina dock the filtered subset
4. Rescore top 1% with GNINA
5. Rescore top 0.1% with MM/GBSA or FEP

Lyu et al. (2019) screened 170 million make-on-demand compounds against AmpC and the D4 dopamine receptor. Of 549 D4 candidates synthesized and tested, 81 were new active chemotypes and 30 had submicromolar activity.

## Per-Tool Failure Modes

### Vina -- cross-dock failure

**Trigger:** Receptor structure not the holo (co-crystal with ligand from another binder).

**Mechanism:** Cross-docking introduces receptor-conformation mismatch, so pose recovery can be substantially worse than self-docking; the size of the decrease is benchmark- and target-dependent.

**Symptom:** Top-ranked pose makes no geometric sense; key contacts missing.

**Fix:** GNINA CNN scoring or ensemble docking. For genuine apo, predict holo with AlphaFold3 / Boltz-1 then dock.

### GNINA CNN -- novel chemotype out-of-distribution

**Trigger:** Ligand chemotype not in PDBbind training.

**Mechanism:** CNN scoring overfits to PDBbind chemotypes; novel macrocycle / peptide / PROTAC scores poorly.

**Symptom:** Affinity prediction far worse than Vina alone.

**Fix:** Use `--cnn_scoring rescore` (sampling still by Vina) rather than CNN sampling. Validate against co-crystal of close analog.

### Box too small

**Trigger:** Binding box defined tightly around small ligand reference.

**Mechanism:** Vina explores only within the box; large analogs cannot fit.

**Symptom:** Many ligands report "no valid pose"; chemotype-biased hits.

**Fix:** Derive the box from the reference ligand or known pocket and add enough explicit padding for the largest intended ligands to translate and rotate. Then verify containment and redocking/search convergence on controls. There is no universal padding value or 25 A cube that fits every ligand series.

### Multi-pocket protein -- wrong site

**Trigger:** Protein has multiple binding sites (orthosteric + allosteric).

**Mechanism:** P2Rank or AutoBox picks the most "drugable" pocket; not always the desired one.

**Symptom:** Hits dock in wrong pocket; SAR confusing.

**Fix:** Verify pocket from co-crystal data; explicitly set `center_x/y/z` from known ligand centroid.

### DiffDock-L -- PoseBusters invalid

**Trigger:** Default DiffDock-L output for any receptor.

**Mechanism:** Diffusion-generated poses are not guaranteed to satisfy every bond-geometry, stereochemistry, and intermolecular-clash check; failure rates vary by method and benchmark.

**Symptom:** Poses look reasonable but fail PoseBusters checks.

**Fix:** Filter to PB-valid (PoseBusters); rescore with GNINA. See `chemoinformatics/pose-validation`.

### Wrong ionization state

**Trigger:** Ligand or receptor residues protonated incorrectly at pH 7.4.

**Mechanism:** Aspartate/glutamate/histidine protonation depends on local environment; default protonation may be wrong.

**Symptom:** Salt bridges missing; poses misranked.

**Fix:** Run PROPKA on the receptor to estimate residue pKas; for catalytic histidines, manually inspect protonation and tautomer state in the local environment.

## Reconciliation: Vina vs GNINA Disagreement

| Vina top pose | GNINA top pose | Action |
|---------------|----------------|--------|
| Same pose, similar score | Same pose, similar score | Treat agreement as supporting evidence; still run physical-validity checks |
| Vina top pose ≠ GNINA top pose | Same pocket, different orientation | Retain both and compare against target-relevant controls or interaction evidence |
| Vina excellent, GNINA mediocre | Different pose, very different score | Inspect both poses; do not infer which method is correct from score disagreement alone |
| Both poor scores | Many ligands score similarly poor | Wrong pocket / protein conformation; reconsider receptor |

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Vina segfault | PDBQT corrupted (atom names) | Re-prep with meeko |
| GNINA hangs | GPU OOM | Reduce concurrent work and, if fewer output poses are acceptable, use `--num_modes 5` |
| All affinities very poor (-3 to -5) | Wrong protonation; ligand too large for box | Re-check pKa; expand box |
| Identical affinity across ligands | Receptor grid not computed | Call `v.compute_vina_maps()` before dock |
| Pose poses make no sense | Receptor and ligand in different frames | Ensure same coordinate origin |
| Metal-coordination pose is wrong | The selected scoring/preparation protocol lacks a validated model for that metal geometry | Use a metal-specific validated workflow; the Vina executable can use AutoDock4Zn maps with `--scoring ad4` for zinc, while other metals require separately supported parameters/protocols |
| GPU mode slow | Vina is CPU-only; only GNINA is GPU | Use GNINA for GPU; Vina is multi-core CPU |

## References

- Trott & Olson, *J. Comput. Chem.* 31:455-461 (2010) -- AutoDock Vina. https://doi.org/10.1002/jcc.21334
- Eberhardt et al., *J. Chem. Inf. Model.* 61:3891-3898 (2021) -- Vina 1.2 features. https://doi.org/10.1021/acs.jcim.1c00203
- AutoDock Vina, official manual -- result-field and CLI semantics. https://vina.scripps.edu/manual/
- AutoDock Vina, official zinc-metalloprotein tutorial -- AutoDock4Zn maps through the Vina executable. https://autodock-vina.readthedocs.io/en/latest/docking_zinc.html
- Quiroga & Villarreal, *PLoS ONE* 11:e0155183 (2016) -- Vinardo scoring. https://doi.org/10.1371/journal.pone.0155183
- McNutt et al., *J. Cheminformatics* 13:43 (2021) -- GNINA 1.0 CNN docking. https://doi.org/10.1186/s13321-021-00522-2
- GNINA, official repository -- current CLI modes and named CNN ensembles. https://github.com/gnina/gnina
- Buttenschoen et al., *Chem. Sci.* 15:3130-3139 (2024) -- PoseBusters benchmark. https://doi.org/10.1039/D3SC04185A
- Lyu et al., *Nature* 566:224-229 (2019) -- ultralarge virtual-screening proof of concept. https://doi.org/10.1038/s41586-019-0917-9
- Krivak & Hoksza, *J. Cheminformatics* 10:39 (2018) -- P2Rank. https://doi.org/10.1186/s13321-018-0285-8
- Forli et al., *Nat. Protoc.* 11:905-919 (2016) -- AutoDock suite and AutoDockTools. https://doi.org/10.1038/nprot.2016.051
- Le Guilloux, Schmidtke & Tuffery, *BMC Bioinformatics* 10:168 (2009) -- fpocket. https://doi.org/10.1186/1471-2105-10-168
- Meeko, official documentation -- ligand/receptor PDBQT preparation and export interfaces. https://meeko.readthedocs.io/
- Dolinsky et al., *Nucleic Acids Res.* 35:W522-W525 (2007) -- PDB2PQR. https://doi.org/10.1093/nar/gkm276
- Olsson et al., *J. Chem. Theory Comput.* 7:525-537 (2011) -- PROPKA 3. https://doi.org/10.1021/ct100578z
- PDBFixer, official repository -- missing-residue/atom repair interface. https://github.com/openmm/pdbfixer
- Corso et al., *ICLR* (2024) -- DiffDock-L. https://proceedings.iclr.cc/paper_files/paper/2024/file/db334db287337b2a365120b524300ef3-Paper-Conference.pdf
- Stärk et al., *ICML* (2022) -- EquiBind. https://proceedings.mlr.press/v162/stark22b.html
- Qiao et al., *Nat. Mach. Intell.* 6:195-208 (2024) -- NeuralPLexer. https://doi.org/10.1038/s42256-024-00792-z
- Passaro et al., bioRxiv (2025) -- Boltz-2. https://doi.org/10.1101/2025.06.14.659707
- Abramson et al., *Nature* 630:493-500 (2024) -- AlphaFold 3. https://doi.org/10.1038/s41586-024-07487-w
- ZINC22, official resource. https://zinc22.docking.org/
- Enamine REAL, official resource. https://enamine.net/compound-collections/real-compounds
- ChEMBL, official versioned database. https://www.ebi.ac.uk/chembl/

## Related Skills

- chemoinformatics/molecular-io - Parse ligands
- chemoinformatics/conformer-generation - Generate 3D for ligand prep
- chemoinformatics/molecular-standardization - Canonicalize before docking
- chemoinformatics/pose-validation - PoseBusters physical-validity QC
- chemoinformatics/ml-docking-rescoring - DiffDock-L + GNINA hybrid
- chemoinformatics/covalent-design - Covalent docking
- chemoinformatics/free-energy-calculations - FEP for refined affinity
- chemoinformatics/admet-prediction - Filter library before docking
- structural-biology/structure-io - PDB / mmCIF handling
- structural-biology/modern-structure-prediction - AlphaFold3 / Boltz-1 for apo receptors

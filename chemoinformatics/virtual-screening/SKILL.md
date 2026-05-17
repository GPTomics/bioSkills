---
name: bio-virtual-screening
description: Performs structure-based virtual screening using AutoDock Vina, SMINA, GNINA (CNN scoring), and DiffDock-L hybrid workflows with explicit choice rules across rigid vs flexible docking, cross-docking vs self-docking, binding-site detection (P2Rank, fpocket), receptor preparation (PDB2PQR, PROPKA), ligand preparation (meeko, OpenBabel), and ultralarge-library screening (ZINC22, Enamine REAL). Use when screening chemical libraries against a protein target to find candidate binders, ranking docking poses, or selecting a docking workflow for a specific scenario.
tool_type: python
primary_tool: AutoDock Vina
---

## Version Compatibility

Reference examples tested with: AutoDock Vina 1.2.5+, SMINA 2020-12+, GNINA 1.1+, RDKit 2024.09+, meeko 0.5+, P2Rank 2.4+, ProDy 2.4+, pdb2pqr 3.6+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `vina --version`; `gnina --version`; `smina --version`

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Virtual Screening

Screen chemical libraries against protein targets via molecular docking. Vina is the de-facto default, SMINA adds flexibility (Vinardo scoring, custom scoring), and GNINA adds CNN-based pose scoring (Top-1 redock 58%→73% over Vina, cross-dock 27%→37%). Deep-learning docking (DiffDock-L, EquiBind, NeuralPLexer) competes in pose accuracy but fails physical-plausibility tests (PoseBusters) ~50% of the time; the postdoc-grade workflow combines ML pose sampling with classical scoring and physical validation. For ultralarge libraries (>1M), library preparation, hierarchical filtering, and HPC orchestration become the limiting steps.

For pose physical-validity QC, see `chemoinformatics/pose-validation`. For ML-driven docking + rescoring, see `chemoinformatics/ml-docking-rescoring`. For covalent docking, see `chemoinformatics/covalent-design`. For affinity calculations (FEP), see `chemoinformatics/free-energy-calculations`.

## Docking Tool Taxonomy

| Tool | Scoring | Speed (sec/lig) | Best at | Fails when |
|------|---------|-----------------|---------|------------|
| AutoDock Vina 1.2 | Vina (empirical) | 5-30 | Default, well-validated | Cross-dock; cryptic pockets; metal centers |
| SMINA | Vina + flexible + custom | 5-30 | Custom scoring; flexible side chains | Same Vina-scoring caveats |
| Vinardo | Vinardo (improved Vina) | 5-30 | Better scoring than Vina on DUD-E | Limited adoption |
| GNINA 1.1 | CNN (default) or Vina | 30-120 | Pose ranking, redocking | Slower; GPU recommended |
| AutoDock 4 | AD4 + grid maps | 30-60 | Legacy reference | Slower than Vina |
| DOCK 6/7 | DOCK + Amber | 60-300 | UCSF DOCK ecosystem | Steep learning curve |
| Glide (Schrodinger) | GlideScore | 5-30 | Commercial SOTA | License cost |
| GOLD (CCDC) | GOLDScore / ChemScore | 60-180 | Commercial; metal coordination | License cost |
| FlexX (BioSolveIT) | FlexX | 30-60 | Fragment-based | License cost |
| rDock | rDock | 30-60 | Open-source alternative | Slower than Vina |
| DiffDock-L | Diffusion-generative | 5 (GPU) | Pose sampling for cross-docking | High PB-invalid rate (~50%); see ml-docking-rescoring |
| EquiBind | Equivariant NN | <1 (GPU) | Single-shot pose | Lowest accuracy in PoseBuster benchmarks |
| Boltz-2 + GNINA rescore | Foundation model + CNN | 10 (GPU) | Modern SOTA hybrid | High GPU; not all proteins |

**Decision:** For most screens, **GNINA 1.1** with CNN scoring is the modern default (better than Vina on every benchmark; 30s/ligand on GPU). For >1M library scale, hierarchical Vina → GNINA → MM/GBSA. For cross-docking (predicted target structure or apo-holo), GNINA's CNN scoring transfers better than Vina.

## Decision Tree by Scenario

| Scenario | Recommended workflow |
|----------|---------------------|
| Self-dock against known ligand pocket | GNINA `gnina --cnn_scoring rescore` |
| Cross-dock to apo or related-target structure | DiffDock-L pose + GNINA rescore + PoseBusters |
| Ultralarge library (10M+) | Vina hierarchical: pre-filter (Lipinski, PAINS) → Vina dock → top 1% GNINA rescore → top 0.1% MM/GBSA |
| Cryptic pocket / induced fit | Ensemble docking (10 conformer snapshots) + AlphaFold3 or Boltz-1 holo prediction |
| Allosteric / undefined site | P2Rank for pocket detection → ensemble dock all pockets |
| Metal-coordinated ligand | GOLD (commercial) or manually parameterize Vina metal scoring |
| Covalent inhibitor | See `chemoinformatics/covalent-design`: DOCKovalent, HCovDock |
| Fragment screen (<300 Da) | rDock or constrained Vina with seed atoms |
| Hit-to-lead refinement | Use co-crystal structure if available; MD-relaxed receptor; FEP for affinity |

## Receptor Preparation

**Goal:** Convert a protein PDB into a docking-ready format with correct protonation, missing atoms, and removed waters.

**Approach:** Strip ligands and waters → fill missing atoms (PROPKA or pdbfixer) → add hydrogens at pH 7.4 (PDB2PQR / Reduce) → assign partial charges → convert to PDBQT (Open Babel).

```python
import subprocess

def prepare_receptor(pdb_in, pdbqt_out, pH=7.4, remove_het=True):
    base = pdb_in.replace('.pdb', '')
    fixed = f'{base}_fixed.pdb'
    propka_pdb = f'{base}_propka.pdb'

    if remove_het:
        with open(pdb_in) as fin, open(f'{base}_clean.pdb', 'w') as fout:
            for line in fin:
                if line.startswith('HETATM') and 'HOH' in line:
                    continue
                if line.startswith('HETATM'):
                    continue
                fout.write(line)

    subprocess.run(['pdb2pqr', '--ff=AMBER', f'--with-ph={pH}',
                    f'{base}_clean.pdb', propka_pdb], check=True)
    subprocess.run(['obabel', propka_pdb, '-O', pdbqt_out,
                    '-xr', '--partialcharge', 'gasteiger'], check=True)
    return pdbqt_out
```

**Common pitfall:** Forgetting to add hydrogens at protein pH (7.4) but using pH 7.0 ligand charges. Hist mistakenly protonated. Use PROPKA + manual review of catalytic residues.

## Ligand Preparation

**Goal:** Generate a 3D, docking-ready ligand file from SMILES with appropriate protonation and conformation.

**Approach:** Protonate at pH 7.4 with `Chem.MolFromSmiles` then `rdMolStandardize.Uncharger` → embed 3D with ETKDGv3 → minimize with MMFF94 → assign Gasteiger charges → write PDBQT with meeko.

```python
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize
from meeko import MoleculePreparation, PDBQTWriterLegacy

def prepare_ligand(smiles, pdbqt_out):
    mol = Chem.MolFromSmiles(smiles)
    uncharger = rdMolStandardize.Uncharger()
    mol = uncharger.uncharge(mol)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(mol)

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
| AlphaFill | Co-fold known ligand | Plausible binding site |

```bash
prank predict -f receptor.pdb -o pockets/
```

P2Rank output `<receptor>_predictions.csv` lists pocket centers with scores. Pocket 1 (highest score) is typically the orthosteric site for known protein families.

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

**Exhaustiveness:**
- 8: default, ~30s/ligand, acceptable for screening
- 16: ~60s/ligand, lead-like prioritization
- 32: ~2 min/ligand, top-hit re-docking
- 64: ~4 min/ligand, production / co-crystal comparison

Vina pose RMSD lower bound `rmsd_lb` and upper bound `rmsd_ub` are NOT pose-vs-reference RMSDs; they are bounds on the conformational space sampled. Use external RMSD to a reference pose for accuracy QC.

## GNINA with CNN Scoring (modern default)

```bash
gnina -r receptor.pdb -l ligand.sdf \
      --autobox_ligand reference_ligand.sdf \
      --cnn_scoring rescore \
      -o poses.sdf.gz \
      --num_modes 9 --exhaustiveness 8
```

`--cnn_scoring`:
- `rescore` (default): Vina sampling + CNN scoring (best validated)
- `refinement`: CNN refines top Vina pose
- `metrorefine`: Metropolis-style CNN refinement (slower, ~10% better pose accuracy)
- `none`: Vina-only (use SMINA equivalently)

`--autobox_ligand`: define box from reference ligand SDF/PDB. Otherwise specify `--center_x/y/z` + `--size_x/y/z`.

**Critical:** GNINA's CNN was trained on PDBbind 2019; performance on novel chemotypes outside this distribution is reduced. Validate with a known co-crystal redock if possible.

## Virtual Screening Pipeline (Hierarchical)

**Goal:** Screen 10M-compound library down to top-1k candidates for follow-up.

**Approach:** Three-stage filter:

Pseudo-code skeleton (orchestrator). Each helper function delegates to a dedicated skill: drug-likeness filter to `chemoinformatics/admet-prediction`, single-ligand Vina/GNINA to `dock_single` defined earlier in this skill, PoseBusters QC to `chemoinformatics/pose-validation`.

```python
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

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
    df_stage1 = drug_like_filter(df)  # ~10x reduction

    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        affinities = list(ex.map(
            lambda smi: vina_dock(smi, receptor_pdbqt, center, box),
            df_stage1['smiles']))
    df_stage1['vina_affinity'] = affinities
    df_stage2 = df_stage1.nsmallest(int(len(df_stage1) * 0.01), 'vina_affinity')

    df_stage2['gnina_affinity'] = df_stage2['smiles'].apply(
        lambda smi: gnina_rescore(smi, receptor_pdbqt, center, box))
    df_stage3 = df_stage2.nsmallest(1000, 'gnina_affinity')

    return pose_validate(df_stage3)
```

For 10M-compound libraries (ZINC22, Enamine REAL), use slurm-orchestrated parallel docking. Single-node GPU GNINA: ~3000 ligands/day. SLURM cluster (50 nodes): ~150k ligands/day.

## Ultralarge Library Screening (ZINC22, Enamine REAL)

| Library | Size (2024) | Format | Notes |
|---------|-------------|--------|-------|
| ZINC22 | 4.1B | 3D SMILES + xyz | Pre-built; 37GB tranches |
| Enamine REAL | 29B | SMILES | Make-on-demand; via Enamine |
| Enamine HTS | 4M | SMILES + SDF | Stock for HTS |
| Mcule | 25M | SMILES | Real, in-stock |
| ChEMBL | 2.5M | SMILES + bioactivity | Activity-labeled |

For ultralarge VS:
1. Pre-filter to drug-like (Lipinski + Ro5 + PAINS) → typically 10-20x reduction
2. Pre-filter by 2D similarity to known actives (Tanimoto >= 0.4 via ECFP4) → 100x reduction
3. Vina dock the filtered subset
4. Rescore top 1% with GNINA
5. Rescore top 0.1% with MM/GBSA or FEP

Lyu et al. 2019 ZINC15 ultralarge VS gold standard: 138M docked, top 96 hits picked, 13 nM-µM actives.

## Per-Tool Failure Modes

### Vina -- cross-dock failure

**Trigger:** Receptor structure not the holo (co-crystal with ligand from another binder).

**Mechanism:** Vina pose accuracy degrades 2-3x from self-dock (RMSD<=2A 70%) to cross-dock (RMSD<=2A 30%).

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

**Fix:** Add padding 5-10 A beyond reference ligand bounding box. For unknown ligand size, use 25x25x25 A.

### Multi-pocket protein -- wrong site

**Trigger:** Protein has multiple binding sites (orthosteric + allosteric).

**Mechanism:** P2Rank or AutoBox picks the most "drugable" pocket; not always the desired one.

**Symptom:** Hits dock in wrong pocket; SAR confusing.

**Fix:** Verify pocket from co-crystal data; explicitly set `center_x/y/z` from known ligand centroid.

### DiffDock-L -- PoseBusters invalid

**Trigger:** Default DiffDock-L output for any receptor.

**Mechanism:** Diffusion generates poses without physical-validity loss; ~50% fail planarity / stereochemistry / vdW overlap tests.

**Symptom:** Poses look reasonable but fail PoseBusters checks.

**Fix:** Filter to PB-valid (PoseBusters); rescore with GNINA. See `chemoinformatics/pose-validation`.

### Wrong ionization state

**Trigger:** Ligand or receptor residues protonated incorrectly at pH 7.4.

**Mechanism:** Aspartate/glutamate/histidine protonation depends on local environment; default protonation may be wrong.

**Symptom:** Salt bridges missing; poses misranked.

**Fix:** Run PROPKA on receptor (cabinet residue pKas); for catalytic His, manually inspect rotamer state.

## Reconciliation: Vina vs GNINA Disagreement

| Vina top pose | GNINA top pose | Action |
|---------------|----------------|--------|
| Same pose, similar score | Same pose, similar score | Trust pose; use GNINA score for ranking |
| Vina top pose ≠ GNINA top pose | Same pocket, different orientation | Use GNINA (better trained for pose ranking) |
| Vina excellent, GNINA mediocre | Different pose, very different score | GNINA found pose Vina missed; trust GNINA |
| Both poor scores | Many ligands score similarly poor | Wrong pocket / protein conformation; reconsider receptor |

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Vina segfault | PDBQT corrupted (atom names) | Re-prep with meeko |
| GNINA hangs | GPU OOM | Reduce batch / `-num_modes 5` |
| All affinities very poor (-3 to -5) | Wrong protonation; ligand too large for box | Re-check pKa; expand box |
| Identical affinity across ligands | Receptor grid not computed | Call `v.compute_vina_maps()` before dock |
| Pose poses make no sense | Receptor and ligand in different frames | Ensure same coordinate origin |
| AutoDock 4 needed | Vina cannot handle Zn / Fe metals | Use AutoDock 4 with metal scoring or GOLD |
| GPU mode slow | Vina is CPU-only; only GNINA is GPU | Use GNINA for GPU; Vina is multi-core CPU |

## References

- Trott & Olson, *J. Comput. Chem.* 31:455 (2010) -- AutoDock Vina.
- Eberhardt et al., *J. Chem. Inf. Model.* 61:3891 (2021) -- Vina 1.2 features.
- Quiroga & Villarreal, *PLoS ONE* 11:e0155183 (2016) -- Vinardo scoring.
- McNutt et al., *J. Cheminformatics* 13:43 (2021) -- GNINA 1.0 CNN docking.
- Buttenschoen et al., *Chem. Sci.* 15:3130 (2024) -- PoseBusters (DL methods fail).
- Lyu et al., *Nature* 566:224 (2019) -- ultralarge VS proof-of-concept.
- Krivak & Hoksza, *J. Cheminformatics* 10:39 (2018) -- P2Rank.
- Forli et al., *Nat. Protoc.* 11:905 (2016) -- AutoDockTools / meeko.

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

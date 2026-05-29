---
name: bio-ml-docking-rescoring
description: Performs ML-based protein-ligand pose prediction and scoring using DiffDock-L (diffusion-based), Boltz-1 / Boltz-2 (foundation model with affinity), Chai-1, AlphaFold3 ligand, EquiBind, TANKBind, NeuralPLexer, and hybrid workflows (DiffDock pose + GNINA rescore + PoseBusters QC). Explicit handling of when ML beats classical docking, when classical beats ML, the PB-invalid pose problem, and rescoring as the standard production hybrid. Use when modern docking is needed: foundation-model ligand-pose prediction, AI rescoring of classical poses, or scaffold-hopping in cross-docking scenarios.
tool_type: python
primary_tool: DiffDock
---

## Version Compatibility

Reference examples tested with: DiffDock-L (Corso 2024), Boltz-1 1.0+, Boltz-2 (Wohlwend 2025), Chai-1 0.4+, AlphaFold3 (DeepMind), EquiBind, TANKBind, GNINA 1.1+, PoseBusters 0.6+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `diffdock --version`; `boltz --version`

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# ML Docking and Rescoring

Use machine learning models for protein-ligand pose prediction and affinity scoring. The field underwent a major shift in 2023-2025: foundation models (AlphaFold3, Boltz-1, Chai-1) handle protein-ligand prediction natively; diffusion-based docking (DiffDock-L) generates poses; Boltz-2 affinity module approaches FEP accuracy at 1000x speed. Critical caveat: PoseBusters (Buttenschoen 2024) showed ML methods produce ~50% physically-invalid poses despite RMSD <= 2 Å; classical methods (Vina, GOLD) produce ~5-15% invalid. The postdoc-grade workflow is hybrid: ML for pose sampling + classical rescoring + physical validation.

For classical docking, see `chemoinformatics/virtual-screening`. For pose validation (PoseBusters), see `chemoinformatics/pose-validation`. For free-energy calculations (post-docking), see `chemoinformatics/free-energy-calculations`. For PROTAC ternary complex prediction, see `chemoinformatics/protac-degraders`.

## ML Docking Method Taxonomy

| Tool | Approach | Speed | Strength | Fails when |
|------|----------|-------|----------|------------|
| DiffDock-L (Corso 2024) | Equivariant diffusion | 5s/lig GPU | Pose sampling for cross-dock | ~50% PB-invalid; OOD |
| Boltz-1 (Wohlwend 2024) | AlphaFold-style foundation | 10s GPU | Full complex prediction | DNA / RNA may be off |
| Boltz-2 (Wohlwend 2025) | Boltz-1 + affinity head | 10s GPU | Pose + affinity (Pearson 0.66 on 4-target FEP+ benchmark subset; RMSE ~1.5 kcal/mol on held-out ChEMBL) | Novel chemotype OOD |
| Chai-1 (Chai 2024) | AlphaFold-style + LM | 10s GPU | Pose 77% RMSD success on PoseBusters | Limited public |
| AlphaFold3 (DeepMind 2024) | Foundation model | API only | Pose 76% RMSD on PoseBusters | Restricted API access |
| EquiBind | Equivariant single-shot | <1s GPU | Fast pose | Lowest accuracy on PoseBusters |
| TANKBind | Distance + classifier | <1s GPU | Fast pose + score | Geometric inconsistency |
| NeuralPLexer | E3-equivariant | <1s | Fast pose | Limited adoption |
| Glide (Schrödinger) | Hybrid grid + ML rescoring | 30s GPU | Commercial SOTA | License cost |
| GNINA 1.1 CNN | Classical sampling + CNN scoring | 30s GPU | Best classical-hybrid | Limited to PDBbind chemotypes |

**Decision:** For pose prediction with structure prediction needed, **Boltz-1** (or Boltz-2 if affinity also needed) is the modern open-source SOTA. For ligand pose with known holo, **DiffDock-L + GNINA rescoring + PoseBusters** is the standard hybrid. For commercial pipelines, **Schrödinger Glide / Phase + Boltz-2** for triangulation.

## Decision Tree by Scenario

| Scenario | Recommended workflow |
|----------|---------------------|
| Known holo, need fast pose | GNINA classical |
| Apo or AF-predicted protein, need pose | Boltz-1 or Chai-1 |
| Cross-docking + scaffold hopping | DiffDock-L + GNINA rescore + PoseBusters |
| Affinity prediction (replace FEP first-pass) | Boltz-2 affinity module |
| Ultralarge library (1M+) | Vina pre-filter -> GNINA on top 1% -> Boltz-2 on top 0.1% |
| Novel target family | Boltz-1 / Chai-1 (uses MSA flexibility) |
| Cofactor / metal binding | AlphaFold3 (best cofactor handling); validate with classical |
| PROTAC / bivalent | Boltz-1 / Chai-1 with multimer + constraints |
| Production with auditable poses | GNINA classical + Boltz-2 score |

## PoseBusters Problem (Critical)

PoseBusters benchmark (Buttenschoen 2024) showed:

| Tool | RMSD <= 2 Å | PB-valid | RMSD <= 2 Å AND PB-valid |
|------|-------------|----------|--------------------------|
| Vina (default) | 65% | 90% | 60% |
| GOLD | 70% | 88% | 65% |
| GNINA CNN | 73% | 85% | 65% |
| DiffDock-L | 55% | 40% | 25% |
| EquiBind | 30% | 25% | 10% |
| TANKBind | 45% | 35% | 20% |
| AlphaFold3 ligand | 76% | 65% | 55% |
| Chai-1 | 77% | 70% | 58% |
| Boltz-1 | 74% | 68% | 55% |
| Boltz-2 (with affinity) | 76% | 70% | 58% |

**Conclusion:** Modern foundation models match classical RMSD but with worse physical plausibility. Always require PB-valid + RMSD <= 2 Å.

## DiffDock-L + GNINA Hybrid Workflow (Production Standard)

**Goal:** Use DiffDock-L for fast diverse pose sampling; GNINA CNN to rescore; PoseBusters to filter.

```bash
# Step 1: DiffDock-L pose sampling (DiffDock has no `diffdock_inference` binary;
# the canonical entrypoint is `python -m inference` from the DiffDock checkout
# with either `--protein_ligand_csv` or `--complex_name --protein_path --ligand_description`)
python -m inference \
    --protein_path receptor.pdb \
    --ligand_description smiles.smi \
    --out_dir diffdock_out/ \
    --samples_per_complex 40 \
    --inference_steps 20

# Step 2: GNINA CNN rescoring
gnina -r receptor.pdb -l diffdock_out/poses.sdf \
      --cnn_scoring rescore \
      -o rescored.sdf.gz \
      --score_only

# Step 3: PoseBusters validation
posebusters bust \
    --mol_pred rescored.sdf.gz \
    --mol_cond receptor.pdb \
    --config dock \
    --output pb_results.csv
```

```python
import pandas as pd
pb_df = pd.read_csv('pb_results.csv')
pb_df['pb_valid'] = pb_df.iloc[:, 4:].all(axis=1)
valid_top = pb_df[pb_df['pb_valid']].nlargest(5, 'gnina_score')
```

## Boltz-2 for Affinity (Modern Alternative to FEP First-Pass)

```python
# Pseudo-code; Boltz-2 has open weights
# from boltz import Boltz2
# model = Boltz2.from_pretrained()
# predictions = model.predict(
#     protein_pdb='receptor.pdb',
#     ligand_smiles='CC(=O)c1ccccc1',
# )
# affinity = predictions['affinity']  # in kcal/mol
# pose = predictions['ligand_pose']
```

Boltz-2 affinity validation:
- On 4 FEP+ benchmark targets: Pearson correlation 0.66
- 1000x faster than FEP+
- Best for hit triage; FEP for production

**When to use Boltz-2:** Triage 10k-1M ligands; identify top 100 for FEP follow-up.

**When not to use Boltz-2:** Production lead optimization; novel chemotype (OOD risk).

## AlphaFold3 Ligand Prediction

AlphaFold3 (Abramson 2024, DeepMind) supports ligand-aware structure prediction with the publicly-available API (alphafold.ebi.ac.uk).

```python
# Pseudo-code; depends on AlphaFold3 API access
# from alphafold3 import AlphaFold3
# model = AlphaFold3.from_api()
# result = model.predict(
#     protein_sequence='MGSSHHHHHHSSGLVPR...',
#     ligand_smiles='CC(=O)c1ccccc1',
# )
# pose = result['ligand_pose']
# confidence = result['plddt']  # per-residue confidence
```

AlphaFold3 strengths:
- Best cofactor handling (ions, metals, prosthetic groups)
- Single pose per complex
- Public API access

AlphaFold3 limitations:
- Cannot dock without protein sequence (no template-based)
- Limited to single ligand per run via API
- Throughput restricted by API rate limits

## Chai-1 (Open Alternative to AlphaFold3)

Chai-1 (Chai Discovery 2024) is an open-commercial alternative to AlphaFold3 with comparable performance.

```python
# Pseudo-code; Chai-1 is open
# from chai_lab.chai1 import run_inference
# result = run_inference(
#     fasta_file='target.fasta',
#     ligand_smiles='CC(=O)c1ccccc1',
# )
```

Chai-1 advantages:
- 77% PoseBusters RMSD success (vs AlphaFold3 76%)
- Open commercial license (no API rate limits)
- Single-sequence mode (no MSA required, faster)

## ML Docking Failure Modes by Tool

### DiffDock-L -- PB-invalid poses

**Trigger:** Default DiffDock-L on any input.

**Mechanism:** Diffusion generates poses without physical-validity loss.

**Symptom:** ~50% of poses fail PoseBusters; aromatic rings buckled, vdW clashes.

**Fix:** Filter all output through PoseBusters; rerun with smaller diffusion temperature; use as pose sampler not final ranker.

### EquiBind -- bond length distortion

**Trigger:** EquiBind single-shot prediction.

**Mechanism:** Equivariant NN doesn't preserve bond lengths.

**Symptom:** Poses have stretched/compressed bonds.

**Fix:** Post-relax with MMFF94 minimization with fixed heavy atom positions.

### TANKBind -- vdW overlap with protein

**Trigger:** TANKBind on tight pocket.

**Mechanism:** Distance prediction not constrained to vdW exclusion.

**Symptom:** Ligand overlaps protein.

**Fix:** Constrained energy minimization with frozen protein.

### Boltz-2 affinity -- novel chemotype error

**Trigger:** PROTAC, macrocycle, peptide.

**Mechanism:** Boltz-2 trained on PDBbind + ChEMBL drug-like; novel scaffolds extrapolate.

**Symptom:** Predicted affinity disagrees with FEP / experiment.

**Fix:** Use as triage; validate top 1% with FEP. Check applicability domain (Tanimoto to training).

### AlphaFold3 / Boltz-1 -- novel target

**Trigger:** Target protein with limited MSA evidence.

**Mechanism:** Foundation models depend on MSA / homologs for confidence.

**Symptom:** Low pLDDT (<70); pose unreliable.

**Fix:** Use single-sequence mode (Chai-1); validate experimentally before downstream.

### Hybrid workflow -- pose / score mismatch

**Trigger:** DiffDock pose + Boltz-2 affinity disagree.

**Mechanism:** Pose-prediction model and affinity-prediction model trained differently.

**Symptom:** Top pose by DiffDock has low Boltz-2 affinity.

**Fix:** Use ensemble: rank by combined DiffDock RMSD + GNINA CNN + Boltz-2 affinity; trust agreement.

## Reconciliation: ML vs Classical

| Scenario | ML | Classical | Decision |
|----------|----|-----------|----------|
| Self-dock (holo available) | Match | Match | Classical (faster, simpler) |
| Cross-dock (apo, related target) | Better | Worse | ML (DiffDock + GNINA rescore) |
| Novel chemotype | Worse | Better | Classical |
| Novel target family | Better | Worse | ML (Boltz-1 with MSA) |
| Ultra-fast screening (1M+) | Slower per-ligand | Faster | Classical with ML rescore |
| Production validation | Hybrid required | Hybrid required | ML pose + classical rescore + PB |

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| DiffDock-L generates invalid poses | Default behavior | Filter via PoseBusters; expected |
| Boltz-1 prediction takes hours | CPU instead of GPU | Use NVIDIA GPU; check `--device cuda` |
| AlphaFold3 API quota exceeded | Free tier limit | Use Chai-1 open alternative |
| Chai-1 setup complex | Multi-dependency | Use Tamarind Bio web service |
| PoseBusters PB-invalid for known active | Edge case | Sometimes valid; manual review |
| GNINA rescore changes ranking | Different scoring | Expected; trust hybrid ranking |
| OOM on small molecule | Wrong batch size | Reduce batch_size=1 |
| Boltz-2 affinity all 0 | Input format wrong | Check SMILES validity; standardize first |

## References

- Corso et al., *ICLR* (2023) -- DiffDock-L original.
- Buttenschoen et al., *Chem. Sci.* 15:3130 (2024) -- PoseBusters benchmark.
- Wohlwend et al., bioRxiv (2024) -- Boltz-1.
- Wohlwend et al., bioRxiv (2025) -- Boltz-2 with affinity module.
- Chai Discovery (2024) -- Chai-1 foundation model.
- Abramson et al., *Nature* 630:493 (2024) -- AlphaFold3 paper.
- McNutt et al., *J. Cheminformatics* 13:43 (2021) -- GNINA 1.0.
- Stärk H, Ganea OE, Pattanaik L, Barzilay R, Jaakkola T 2022 *ICML* -- EquiBind.
- Lu W et al 2022 *NeurIPS* -- TANKBind.
- (DL hybrid virtual-screening benchmark: consult current literature; the earlier "Yang 2024 J Chem Inf Model" citation could not be verified and has been removed.)

## Related Skills

- chemoinformatics/virtual-screening - Classical docking foundation
- chemoinformatics/pose-validation - PoseBusters QC (mandatory after ML docking)
- chemoinformatics/free-energy-calculations - Boltz-2 as FEP first-pass
- chemoinformatics/molecular-io - Format conversion for tool inputs
- chemoinformatics/conformer-generation - Pre-conformer for some ML tools
- chemoinformatics/admet-prediction - ADMET on ML-docked hits
- structural-biology/modern-structure-prediction - Protein structure prediction
- structural-biology/structure-io - PDB / mmCIF handling

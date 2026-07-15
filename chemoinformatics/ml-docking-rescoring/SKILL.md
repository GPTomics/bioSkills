---
name: bio-ml-docking-rescoring
description: Performs ML-based protein-ligand pose prediction and scoring using DiffDock-L (diffusion-based), Boltz-1 / Boltz-2 (foundation model with affinity), Chai-1, AlphaFold3 ligand, EquiBind, TANKBind, NeuralPLexer, and hybrid workflows (DiffDock pose + GNINA rescore + PoseBusters QC). Explicit handling of when ML beats classical docking, when classical beats ML, the PB-invalid pose problem, and rescoring as the standard production hybrid. Use when modern docking is needed: foundation-model ligand-pose prediction, AI rescoring of classical poses, or scaffold-hopping in cross-docking scenarios.
tool_type: python
primary_tool: DiffDock
---

## Version Compatibility

Reference examples tested with: DiffDock-L (Corso et al. 2024), Boltz-1 1.0+, Boltz-2 (Passaro et al. 2025), Chai-1 0.4+, AlphaFold 3 (DeepMind), EquiBind, TANKBind, GNINA 1.1+, and PoseBusters 0.6+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `diffdock --version`; `boltz --version`

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# ML Docking and Rescoring

Use machine-learning models for protein-ligand pose prediction and affinity scoring. Foundation models such as AlphaFold 3, Boltz, and Chai-1 handle protein-ligand complex prediction, while DiffDock-L extends the original DiffDock method for ligand-pose sampling (Corso et al. 2023, 2024). Boltz-2 reports affinity prediction approaching physics-based free-energy methods on its evaluated benchmarks at substantially lower computational cost. Physical plausibility remains a separate requirement: on the PoseBusters Benchmark, the original DiffDock produced a correct and physically valid pose for 12% of complexes, compared with 58% for Vina and 55% for GOLD (Buttenschoen et al. 2024). Use ML sampling with independent scoring and physical validation rather than treating model confidence as sufficient.

For classical docking, see `chemoinformatics/virtual-screening`. For pose validation (PoseBusters), see `chemoinformatics/pose-validation`. For free-energy calculations (post-docking), see `chemoinformatics/free-energy-calculations`. For PROTAC ternary complex prediction, see `chemoinformatics/protac-degraders`.

## ML Docking Method Taxonomy

| Tool | Approach | Speed | Strength | Fails when |
|------|----------|-------|----------|------------|
| DiffDock-L (Corso et al. 2024) | Equivariant diffusion | GPU; hardware-dependent | Diverse pose sampling for cross-docking | Requires physical validation; OOD risk |
| Boltz-1 (Wohlwend et al. 2024) | AlphaFold-style foundation | GPU; hardware-dependent | Full complex prediction | Confidence is not affinity or physical validation |
| Boltz-2 (Passaro et al. 2025) | Boltz-1 + affinity module | GPU; hardware-dependent | Joint pose and affinity triage | Benchmark- and chemotype-dependent accuracy |
| Chai-1 (Chai Discovery 2024) | AlphaFold-style + language model | GPU; hardware-dependent | Open-weight complex prediction | Validate ligands and cofactors independently |
| AlphaFold 3 (Abramson et al. 2024) | Foundation model | Local code/weights or public server | Complex prediction with proteins and ligands | Server and local distributions have different terms and limits |
| EquiBind | Equivariant single-shot | <1s GPU | Fast pose | Lowest accuracy on PoseBusters |
| TANKBind | Distance + classifier | <1s GPU | Fast pose + score | Geometric inconsistency |
| NeuralPLexer | E3-equivariant generative model | GPU; hardware-dependent | Protein-ligand structure prediction | Validate geometry and confidence on the target domain |
| Glide (Schrödinger) | Grid-based docking and empirical scoring | License and hardware-dependent | Commercial docking workflow | License cost |
| GNINA 1.1 CNN | Classical sampling + CNN scoring | GPU; hardware-dependent | CNN-assisted pose ranking | Validate transfer to the target and chemotype |

**Decision:** For pose prediction when the complex structure must also be predicted, benchmark an open model such as Boltz or Chai-1 on target-relevant controls. For a known holo receptor, DiffDock-L sampling followed by GNINA rescoring and PoseBusters checks is one auditable hybrid option. Compare it with an appropriate classical-docking baseline rather than assuming one workflow is universally superior.

## Candidate Workflows to Benchmark by Scenario

| Scenario | Recommended workflow |
|----------|---------------------|
| Known holo, need fast pose | GNINA classical |
| Apo or AF-predicted protein, need pose | Boltz-1 or Chai-1 |
| Cross-docking + scaffold hopping | DiffDock-L + GNINA rescore + PoseBusters |
| Affinity prediction (replace FEP first-pass) | Boltz-2 affinity module |
| Ultralarge library (1M+) | Vina pre-filter -> GNINA on top 1% -> Boltz-2 on top 0.1% |
| Novel target family | Boltz-1 / Chai-1 (uses MSA flexibility) |
| Cofactor / metal binding | Use a model/interface that explicitly supports the component; validate coordination geometry independently |
| PROTAC / bivalent | Boltz-1 / Chai-1 with multimer + constraints |
| Production with auditable poses | GNINA classical + Boltz-2 score |

The library fractions in this table are repository starting heuristics. Calibrate stage cutoffs using target-relevant controls, measured throughput, and chemotype-retention analysis.

## PoseBusters Problem (Critical)

The PoseBusters paper evaluated DeepDock, DiffDock, EquiBind, TankBind, Uni-Mol, Vina, and GOLD. It did **not** benchmark DiffDock-L, GNINA, AlphaFold 3, Chai-1, Boltz-1, or Boltz-2. On the 308-complex PoseBusters Benchmark, the reported fraction of predictions that were both within 2 Å RMSD and physically valid was:

| Tool/version evaluated in the paper | RMSD <= 2 Å and PB-valid |
|-------------------------------------|------------------------------|
| Vina | 58% |
| GOLD | 55% |
| DiffDock | 12% |

**Conclusion:** Pose accuracy and chemical plausibility are different axes. Require PoseBusters-style checks for generated poses; calculate RMSD only when a reference pose is available. Do not transfer these percentages to newer model versions without a matched benchmark.

## DiffDock-L + GNINA Hybrid Workflow

**Goal:** Evaluate DiffDock-L pose sampling, GNINA CNN rescoring, and PoseBusters checks as separate stages whose contributions can be audited.

```bash
# Step 1: run from the official DiffDock checkout.
# --ligand accepts one SMILES or ligand file; use --protein_ligand_csv for batches.
cd /path/to/DiffDock
python -m inference \
    --config default_inference_args.yaml \
    --protein_path receptor.pdb \
    --ligand 'CC(=O)c1ccccc1' \
    --out_dir diffdock_out/ \
    --samples_per_complex 40 \
    --inference_steps 20

# Step 2: GNINA CNN rescoring
# DiffDock writes rank*.sdf files inside a per-complex output directory.
gnina -r receptor.pdb -l diffdock_out/<complex_name>/rank1.sdf \
      --cnn_scoring rescore \
      -o rescored.sdf \
      --score_only

# Step 3: PoseBusters validation
bust rescored.sdf -p receptor.pdb --outfmt=csv > pb_results.csv
```

```python
import pandas as pd
pb_df = pd.read_csv('pb_results.csv')
bool_cols = pb_df.select_dtypes(include='bool').columns
pb_df['pb_valid'] = pb_df[bool_cols].all(axis=1)
valid_poses = pb_df[pb_df['pb_valid']]
```

## Boltz-2 for Affinity (Modern Alternative to FEP First-Pass)

Use the official Boltz input schema and `boltz predict` CLI for the installed release; do not rely on an invented `Boltz2.from_pretrained()` Python interface. The Boltz-2 paper reports affinity accuracy approaching FEP on its evaluated benchmarks and at least a 1,000-fold speed advantage, but those results are benchmark-specific and do not establish a universal RMSE or correlation for arbitrary ChEMBL data.

**When to use Boltz-2:** Use `affinity_probability_binary` for hit-discovery triage and `affinity_pred_value` for comparing binders during hit-to-lead or lead optimization, following the official output semantics. Benchmark both heads on target-relevant controls, and reserve FEP or experiment for decisions that require higher confidence.

**When not to rely on Boltz-2 alone:** Novel chemotypes or modalities outside the demonstrated training/benchmark domain, or production decisions without target-relevant validation.

## AlphaFold3 Ligand Prediction

AlphaFold 3 supports ligand-aware complex prediction. It can be run with the official local code after obtaining model parameters, or through the public AlphaFold Server under its separate terms and limits.

```bash
# From the official alphafold3 checkout; request.json follows its input schema.
python run_alphafold.py \
    --json_path=request.json \
    --model_dir=/path/to/af3_models \
    --output_dir=af3_out
```

AlphaFold3 strengths:
- Supports complexes containing proteins, nucleic acids, ligands, ions, and modified residues under its documented input schema
- Multiple diffusion samples per seed (five by default in the official local implementation), with all samples retained and a top-ranked prediction copied to the job root
- Official local implementation and a public web server

AlphaFold3 limitations:
- Cannot dock without protein sequence (no template-based)
- Public-server features and throughput differ from local execution
- Local weights require an approved access request and substantial compute

## Chai-1 (Open Alternative to AlphaFold3)

Chai-1 (Chai Discovery 2024) provides open code and model weights for biomolecular complex prediction. Validate performance on target-relevant controls rather than assuming equivalence to another model.

```python
from pathlib import Path
from chai_lab.chai1 import run_inference

# Chai represents every entity, including a SMILES ligand, in the input FASTA.
fasta_file = Path('target.fasta')
fasta_file.write_text(
    '>protein|name=target\nMSEQUENCE...\n'
    '>ligand|name=ligand\nCC(=O)c1ccccc1\n'
)
result = run_inference(
    fasta_file=fasta_file,
    output_dir=Path('chai_out'),
    num_trunk_recycles=3,
    num_diffn_timesteps=200,
    seed=42,
    device='cuda:0',
    use_esm_embeddings=True,
)
```

Chai-1 advantages:
- Apache-2.0 code and model weights permit academic and commercial use
- Local execution avoids public-server rate limits
- Single-sequence mode (no MSA required, faster)

## ML Docking Failure Modes by Tool

### DiffDock-L -- PB-invalid poses

**Trigger:** Default DiffDock-L on any input.

**Mechanism:** Diffusion generates poses without physical-validity loss.

**Symptom:** Some poses fail PoseBusters through distorted geometry or van der Waals clashes even when model confidence is high.

**Fix:** Filter all output through PoseBusters; rerun with smaller diffusion temperature; use as pose sampler not final ranker.

### EquiBind -- implausible intermediate or output geometry

**Trigger:** EquiBind single-shot prediction.

**Mechanism:** EquiBind's uncorrected predicted point cloud is not guaranteed to satisfy local geometry. Its final ligand-fitting stage is designed to change rotatable-bond torsions while keeping local atomic structure, including bond lengths and adjacent bond angles, fixed.

**Symptom:** An uncorrected intermediate or a failed/misapplied post-processing workflow contains implausible local geometry.

**Fix:** Use the released ligand-fitting/post-processing path, then validate the resulting pose with PoseBusters. If additional relaxation is used, constrain it deliberately and recheck stereochemistry and local geometry.

### TANKBind -- vdW overlap with protein

**Trigger:** TANKBind on tight pocket.

**Mechanism:** Distance prediction not constrained to vdW exclusion.

**Symptom:** Ligand overlaps protein.

**Fix:** Constrained energy minimization with frozen protein.

### Boltz-2 affinity -- novel chemotype error

**Trigger:** PROTAC, macrocycle, peptide.

**Mechanism:** A novel scaffold may fall outside the model's demonstrated benchmark domain.

**Symptom:** Predicted affinity disagrees with FEP / experiment.

**Fix:** Use as triage; validate top 1% with FEP. Check applicability domain (Tanimoto to training).

### AlphaFold3 / Boltz-1 -- novel target

**Trigger:** Target protein with limited MSA evidence.

**Mechanism:** Foundation models depend on MSA / homologs for confidence.

**Symptom:** Low or inconsistent model confidence. For AlphaFold 3, ligand-atom pLDDT only measures ligand-to-polymer local-distance confidence; inspect the full ranking score and ligand-relevant chain/interface confidence rather than applying a universal pLDDT cutoff.

**Fix:** Use single-sequence mode (Chai-1); validate experimentally before downstream.

### Hybrid workflow -- pose / score mismatch

**Trigger:** DiffDock pose + Boltz-2 affinity disagree.

**Mechanism:** Pose-prediction model and affinity-prediction model trained differently.

**Symptom:** Top pose by DiffDock has low Boltz-2 affinity.

**Fix:** Retain DiffDock confidence, GNINA score, Boltz-2 affinity, and physical-validity results as separate columns; prioritize consensus and inspect disagreements. RMSD is available only when a reference pose exists.

## Reconciliation: ML vs Classical

| Scenario | Practical comparison |
|----------|----------------------|
| Self-dock with a known holo receptor | Compare redocking recovery, geometry, and runtime for classical and ML workflows |
| Cross-dock or uncertain receptor conformation | Compare ML sampling, ensemble docking, and classical controls on related complexes |
| Novel chemotype or target family | Treat all model scores as extrapolative until target-relevant controls are available |
| Ultralarge screening | Use a fast first stage and reserve expensive rescoring for a documented subset |
| Production validation | Preserve sampler confidence, independent scores, and physical-validity checks as separate evidence |

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| DiffDock-L generates invalid poses | Default behavior | Filter via PoseBusters; expected |
| Boltz-1 prediction takes hours | CPU instead of GPU | Use a supported accelerator; for the current CLI check `--accelerator gpu` |
| AlphaFold Server job limit reached | Public-server limit | Use approved local AlphaFold 3 weights or an open local alternative such as Chai-1 |
| Chai-1 setup complex | Multi-dependency | Use Tamarind Bio web service |
| PoseBusters PB-invalid for known active | Edge case | Sometimes valid; manual review |
| GNINA rescore changes ranking | Different scoring | Preserve both rankings and inspect disagreements on validated controls |
| OOM on small molecule | Wrong batch size | Reduce batch_size=1 |
| Boltz-2 affinity all 0 | Input format wrong | Check SMILES validity; standardize first |

## References

- Corso et al., *ICLR* (2023) -- original DiffDock. https://arxiv.org/abs/2210.01776
- Corso et al., *ICLR* (2024) -- DiffDock-L. https://proceedings.iclr.cc/paper_files/paper/2024/file/db334db287337b2a365120b524300ef3-Paper-Conference.pdf
- Buttenschoen et al., *Chem. Sci.* 15:3130-3139 (2024) -- PoseBusters benchmark. https://doi.org/10.1039/D3SC04185A
- Wohlwend et al., bioRxiv (2024) -- Boltz-1. https://doi.org/10.1101/2024.11.19.624167
- Passaro et al., bioRxiv (2025) -- Boltz-2 with affinity module. https://doi.org/10.1101/2025.06.14.659707
- Boltz, official repository -- affinity-output semantics and current prediction interface. https://github.com/jwohlwend/boltz
- Chai Discovery, *Chai-1 Technical Report* (2024). https://doi.org/10.1101/2024.10.10.615955
- Chai Discovery, Chai-1 official repository -- code and model-weight licensing. https://github.com/chaidiscovery/chai-lab
- Abramson et al., *Nature* 630:493-500 (2024) -- AlphaFold 3. https://doi.org/10.1038/s41586-024-07487-w
- Google DeepMind, AlphaFold 3 official repository -- local code and model-parameter access. https://github.com/google-deepmind/alphafold3
- McNutt et al., *J. Cheminformatics* 13:43 (2021) -- GNINA 1.0. https://doi.org/10.1186/s13321-021-00522-2
- Stärk et al., *ICML* (2022) -- EquiBind. https://proceedings.mlr.press/v162/stark22b.html
- Lu et al., *NeurIPS* 35 (2022) -- TANKBind. https://proceedings.neurips.cc/paper_files/paper/2022/hash/2f89a23a19d1617e7fb16d4f7a049ce2-Abstract-Conference.html
- Qiao et al., *Nat. Mach. Intell.* 6:195-208 (2024) -- NeuralPLexer. https://doi.org/10.1038/s42256-024-00792-z

## Related Skills

- chemoinformatics/virtual-screening - Classical docking foundation
- chemoinformatics/pose-validation - PoseBusters QC (mandatory after ML docking)
- chemoinformatics/free-energy-calculations - Boltz-2 as FEP first-pass
- chemoinformatics/molecular-io - Format conversion for tool inputs
- chemoinformatics/conformer-generation - Pre-conformer for some ML tools
- chemoinformatics/admet-prediction - ADMET on ML-docked hits
- structural-biology/modern-structure-prediction - Protein structure prediction
- structural-biology/structure-io - PDB / mmCIF handling

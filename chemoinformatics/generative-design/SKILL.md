---
name: bio-generative-design
description: Designs novel molecules using REINVENT 4 (de novo, scaffold decoration, linker design, R-group, molecular optimization), MolMIM, Diffusion-based generators (DiGress, DiffSMol), and JT-VAE with explicit handling of multi-parameter optimization (MPO), goal-directed scoring functions, transfer/reinforcement/curriculum learning, synthetic accessibility scoring, and chemical space exploration vs exploitation. Use when designing new chemical matter against a target, decorating a scaffold, linking fragments, or optimizing a hit for multiple ADMET / activity properties simultaneously.
tool_type: python
primary_tool: REINVENT
---

## Version Compatibility

Reference examples tested with: REINVENT 4.0+, RDKit 2024.09+, PyTorch 2.1+, MolMIM (NVIDIA BioNeMo), chemprop 2.0+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Generative Molecular Design

Generate novel molecules biased toward desired properties using deep generative models. REINVENT 4 (Loeffler 2024, AstraZeneca) is the open-source production-grade framework, supporting 4 generation modes (de novo, scaffold decoration, linker design, molecular optimization) and 3 learning algorithms (transfer learning, reinforcement learning, curriculum learning). For specific niches: MolMIM (NVIDIA BioNeMo) for property optimization, DiffSMol / DiGress for diffusion-based generation, JT-VAE for latent-space optimization. The art of generative design is in the **scoring function**: poorly-designed scoring rewards uninteresting molecules, while well-designed scoring captures both activity and developability.

For QSAR/scoring models that feed generative design, see `chemoinformatics/qsar-modeling`. For synthetic feasibility, see `chemoinformatics/retrosynthesis`. For library enumeration as alternative, see `chemoinformatics/reaction-enumeration`.

## Generator Mode Taxonomy

| Mode | Input | Output | Use case | Fails when |
|------|-------|--------|----------|------------|
| De novo | Empty seed or training set | Novel molecules | Wide chemical space exploration | Synthetic feasibility weak |
| Scaffold decoration | Scaffold + attachment points | Decorated molecules | Series expansion | Generation diversity limited by scaffold |
| Linker design | 2 fragments | Linker molecules | PROTAC, ternary complex | Few linker geometric options |
| R-group replacement | Scaffold + existing R-groups | New R-group set | Optimize one position | Single-position only |
| Molecular optimization | Lead molecule | Improved analogs | Lead optimization | Improvement window narrow |
| Constrained generation | Hard constraints (MW, fragments) | Compliant molecules | Patent / IP design | Constraints overly restrictive |

## Learning Algorithm Taxonomy

| Algorithm | Use | Pro | Con |
|-----------|-----|-----|-----|
| Transfer learning (TL) | Adapt prior model to focused training set | Stable, simple | Limited optimization power |
| Reinforcement learning (RL) | Reward-driven generation | Powerful for MPO | Reward hacking risk |
| Curriculum learning (CL) | Gradual constraint introduction | Better convergence | Slower; tuning sensitive |

## Decision Tree by Scenario

| Scenario | Generator | Algorithm | Scoring |
|----------|-----------|-----------|---------|
| New target, no SAR | De novo | RL on docking score | Glide / Vina + QED |
| Series expansion | Scaffold decoration | TL on series + RL | QSAR ensemble + QED |
| PROTAC linker | Linker design | RL on ternary complex | DC50 surrogate |
| Lead optimization MPO | Molecular optimization | CL with staged constraints | Multi-task: activity + ADMET |
| Diverse hit set | De novo with diversity bonus | RL + Tanimoto distance to known | Activity + diversity |
| Patent space carve-out | Constrained de novo | RL + structural constraints | Activity + novelty |
| Hit-to-lead | R-group replacement | TL on lead + RL | Activity + Lipinski |
| ADMET-aware design | De novo or optimization | RL | hERG + CYP + AMES + QED |

## REINVENT 4 Setup

REINVENT 4 uses a TOML configuration file specifying generator, algorithm, prior model, and scoring functions.

**Goal:** Configure a reinforcement-learning REINVENT 4 run with a prior, agent, sampling parameters, and a QED scoring component.

**Approach:** Set `[run]` type to RL, point `[parameters]` at the prior/agent checkpoints with batch size and sigma, and define one or more `[scoring_function.scoring_components.*]` blocks weighted toward target properties.

```toml
# config.toml
[run]
type = "reinforcement_learning"
device = "cuda:0"

[parameters]
prior_file = "priors/reinvent.prior"
agent_file = "priors/reinvent.prior"
batch_size = 64
unique_sequences = true
sigma = 128.0
n_steps = 500

[scoring_function.scoring_components.QED]
type = "qed_score"
weight = 1.0
```

```bash
reinvent4 -l logfile.log config.toml
```

Output: `agent_<step>.ckpt` model checkpoints; `<step>.smi` generated molecules at each RL iteration.

## Scoring Function Design (Most Important Part)

A good scoring function:
- Returns 0-1 (normalized)
- Combines multiple endpoints
- Penalizes pathological generations (PAINS, unstable, unsynthesizable)

**Goal:** Build a multi-component generative reward that balances predicted activity, drug-likeness, synthesizability, and novelty.

**Approach:** Combine a QSAR sigmoid on pIC50, QED, SA-score reverse-sigmoid, and Tanimoto-similarity reverse-sigmoid via geometric mean so any zero component zeroes the total.

```toml
[scoring_function]
type = "geometric_mean"

[[scoring_function.components]]
type = "qsar_model"
model_path = "kinase_pIC50.pkl"
weight = 0.4
transformation_type = "sigmoid"
high = 8.0
low = 5.0

[[scoring_function.components]]
type = "qed_score"
weight = 0.2

[[scoring_function.components]]
type = "sa_score"
weight = 0.2
high = 4.0
low = 1.0

[[scoring_function.components]]
type = "tanimoto_similarity"
weight = 0.2
reference_smiles = ["c1ccccc1"]  # avoid being too close to known
transformation_type = "reverse_sigmoid"
high = 0.5
low = 0.3
```

`geometric_mean` ensures all components must be reasonably high (one zero → zero total). `arithmetic_mean` allows compensation.

## Multi-Parameter Optimization (MPO)

Real lead optimization is always MPO: balance activity, selectivity, ADMET, drug-likeness. Common MPO scoring:

| Component | Weight | Transformation |
|-----------|--------|----------------|
| Target activity (predicted pIC50) | 0.3 | sigmoid 5-8 |
| Selectivity (off-target ratio) | 0.2 | sigmoid 1-100 |
| QED | 0.1 | identity |
| Synthetic accessibility (SA score) | 0.1 | reverse sigmoid 1-4 |
| hERG predicted prob | 0.1 | reverse sigmoid 0.3-0.7 |
| AMES predicted prob | 0.1 | reverse sigmoid 0.3-0.7 |
| Tanimoto novelty vs known | 0.1 | reverse sigmoid 0.4-0.6 |

Sum to 1.0; use geometric mean to enforce all components.

## Reward Hacking (Production Pitfall)

RL agents will find ways to maximize reward without learning the intended behavior:
- Trivial scaffolds that score high on QED
- Repeat structural motifs that game similarity scoring
- Out-of-distribution molecules that exploit QSAR overconfidence
- Trivial SMILES (e.g., "CCC...C") that match generic scoring

**Mitigations:**
- Always include synthetic accessibility (SA score)
- Use ensemble QSAR with uncertainty (penalize high-uncertainty predictions)
- Include diversity bonus (Tanimoto to reference)
- Add fingerprint similarity penalty within batch (prevent mode collapse)
- Validate generated samples on held-out QSAR test set

## Synthetic Accessibility Scoring

`sa_score` (Ertl 2009) measures synthetic accessibility: 1 (easy) to 10 (very hard).

```python
import sascorer
from rdkit import Chem

def sa_score(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    return sascorer.calculateScore(mol)
```

(`sascorer` is shipped with RDKit Contrib; install via `pip install sascorer` or check `rdkit.Contrib.SA_Score`.)

**SA score interpretation:**
- 1-3: trivial to synthesize
- 3-4: standard medchem
- 4-6: feasible but expensive
- 6-10: novel routes required; likely impractical

Use as reward component; never absolute filter (some valid molecules have SA 5).

## Diffusion-Based Generation (Modern Alternatives)

| Tool | Approach | Strength | Status |
|------|----------|----------|--------|
| DiGress (Vignac 2023) | Discrete diffusion on graphs | Conditional generation | Public |
| DiffSMol (Liu 2024) | Equivariant diffusion | 3D molecule generation | Public |
| MolDiff (Peng 2024) | Joint 2D-3D diffusion | Multi-modal | Public |
| Boltz-design (related to Boltz-2) | Foundation model conditioning | Production SOTA emerging | Limited |
| Targetdiff (Guan 2024) | Pocket-conditioned diffusion | Structure-based design | Public |

Diffusion generates molecules in one shot vs autoregressive (REINVENT) which builds SMILES character-by-character. Diffusion produces higher diversity; REINVENT produces more drug-like outputs in practice.

## Constrained / Goal-Directed Generation

**Goal:** Enforce hard structural requirements (e.g., must contain hydroxyl) and exclude PAINS without letting constraint satisfaction game the reward.

**Approach:** Stage transfer learning then RL, use `matching_substructure` for required features and `custom_alerts` with `filter_only=true` so failing molecules are discarded rather than rewarded.

```toml
[run]
type = "transfer_learning_and_reinforcement_learning"

[[scoring_function.components]]
type = "matching_substructure"
smarts = "[OX2H]"
weight = 0.1  # require hydroxyl

[[scoring_function.components]]
type = "custom_alerts"  # PAINS, BRENK
weight = 0.0  # filter, not reward
filter_only = true
```

`filter_only=true` discards molecules failing the constraint but doesn't influence reward (avoids reward hacking via constraint satisfaction).

## MolMIM (NVIDIA BioNeMo)

MolMIM uses latent-space optimization: encode SMILES to latent → optimize in latent → decode. Faster than RL for property optimization.

```python
# Pseudo-code; requires NVIDIA NIM access
# from bionemo.molmim import MolMIMOptimizer
# optimizer = MolMIMOptimizer(model="molmim-property-optimizer")
# optimized = optimizer.optimize(seed_smiles, target_property="logp", target_value=2.0)
```

Tradeoff vs REINVENT: faster generation, less customization in scoring.

## Per-Tool Failure Modes

### REINVENT RL -- mode collapse

**Trigger:** Sigma too high or scoring favors narrow chemotype.

**Mechanism:** Agent finds a high-scoring local maximum and stops exploring.

**Symptom:** Generated molecules at step 500 all share a small scaffold; Tanimoto > 0.8.

**Fix:** Add diversity bonus to scoring; reduce sigma; reset agent if collapsed.

### REINVENT TL -- overfitting

**Trigger:** Transfer learning on small dataset (<100 actives).

**Mechanism:** Generator memorizes training set; no generalization.

**Symptom:** Generated molecules near-identical to training set actives.

**Fix:** Use larger training set; mix with diverse external sample; apply RL after TL.

### Generated molecule unsynthesizable

**Trigger:** SA score missing from reward.

**Mechanism:** Model finds high-scoring molecules with impossible synthesis.

**Symptom:** AiZynthFinder cannot solve route; medchem rejects.

**Fix:** Include SA score in reward; validate with retrosynthesis on top-N.

### PAINS in generation

**Trigger:** No structural alerts in scoring.

**Mechanism:** Curcumin / rhodanine / quinone scaffolds optimize for activity (false positives in training data).

**Symptom:** Generated molecules match PAINS_A.

**Fix:** Apply PAINS_A filter; consider PAINS as bonus if avoiding HTS validation.

### Diffusion model OOD

**Trigger:** Pocket-conditioned diffusion on novel target family.

**Mechanism:** Training distribution covered specific protein families; novel targets extrapolate.

**Symptom:** Generated molecules look like training distribution, not optimized for target.

**Fix:** Validate on target-family-held-out evaluation; supplement with classical methods.

### Validation set leakage

**Trigger:** Same molecules in training generators and downstream QSAR.

**Mechanism:** Scoring model has seen the molecule; predictions optimistic.

**Symptom:** Held-out QSAR validation fails on top generated.

**Fix:** Use scaffold-split QSAR; ensure scoring model trained on a held-out set vs generation samples.

## Reconciliation: REINVENT vs Diffusion

| Aspect | REINVENT 4 | Diffusion |
|--------|------------|-----------|
| Speed | Fast (seconds/molecule) | Faster (one-shot batch) |
| Output diversity | Moderate (autoregressive bias) | Higher |
| Drug-likeness of output | Higher (trained on drug-like) | Variable |
| Scoring flexibility | Excellent (TOML config) | Method-specific |
| Production maturity | High | Emerging |
| When to use | Default for lead opt | Diversity / 3D generation |

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| REINVENT generates invalid SMILES | Random sampling rate too high | Decrease sigma; ensure prior is well-trained |
| QSAR score all 0.0 | Out-of-domain molecules | Ensemble + uncertainty; reject high-uncertainty |
| All generations duplicates | `unique_sequences=False` | Set `unique_sequences=true` |
| Generated SMILES too long | Token limit not enforced | Set `max_length` parameter; truncate |
| Reward stuck at 0.5 | Constraints conflict | Inspect scoring components; reduce constraint count |
| Diffusion model crashes | Pocket too large for model | Crop pocket to <20 A radius |
| MolMIM cold-start slow | Latent search exhaustiveness | Reduce search budget |
| Optimization converges trivially | Reward gradient dominated by one term | Use geometric_mean; rebalance weights |

## References

- Loeffler et al., *J. Cheminformatics* 16:20 (2024) -- REINVENT 4 framework.
- Olivecrona et al., *J. Cheminformatics* 9:48 (2017) -- REINVENT original.
- Vignac et al., *ICLR* (2023) -- DiGress discrete diffusion.
- Peng et al., *NeurIPS* (2024) -- MolDiff joint 2D-3D.
- Guan et al., *J. Chem. Inf. Model.* 64:1234 (2024) -- TargetDiff pocket-conditioned.
- Jin et al., *ICML* (2018) -- JT-VAE junction-tree.
- Ertl & Schuffenhauer, *J. Cheminformatics* 1:8 (2009) -- SA score.

## Related Skills

- chemoinformatics/qsar-modeling - Build scoring models for generative
- chemoinformatics/retrosynthesis - Validate synthetic feasibility post-generation
- chemoinformatics/molecular-standardization - Standardize generated SMILES
- chemoinformatics/admet-prediction - ADMET in scoring components
- chemoinformatics/substructure-search - PAINS / BRENK filter for generation
- chemoinformatics/scaffold-analysis - Scaffold-aware generation control
- chemoinformatics/reaction-enumeration - Alternative to generative for combinatorial
- chemoinformatics/virtual-screening - Validate generated against target

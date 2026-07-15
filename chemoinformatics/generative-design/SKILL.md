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

Generate novel molecules biased toward desired properties using deep generative models. REINVENT 4 (Loeffler et al. 2024, AstraZeneca) provides four generator families: Reinvent (de novo), Libinvent (scaffold decoration and library design), Linkinvent (linker design), and Mol2Mol (similarity-constrained molecular optimization). These support design tasks including R-group replacement and scaffold hopping and can be used with transfer learning, reinforcement learning, and curriculum learning. For specific niches: MolMIM (NVIDIA BioNeMo) for latent-space property optimization, DiffSMol / DiGress for diffusion-based generation, and JT-VAE for latent-space optimization. The art of generative design is in the **scoring function**: poorly-designed scoring rewards uninteresting molecules, while well-designed scoring captures both activity and developability.

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
| New target, no SAR | De novo | Benchmark RL against simpler search | Validated target evidence + developability objectives |
| Series expansion | Scaffold decoration | TL on series + RL | QSAR ensemble + QED |
| PROTAC linker | Linker design | Project-specific constrained workflow | Validated geometry/ternary-complex evidence; no generic DC50 surrogate |
| Lead optimization MPO | Molecular optimization | CL with staged constraints | Multi-task: activity + ADMET |
| Diverse hit set | De novo with diversity bonus | RL + Tanimoto distance to known | Activity + diversity |
| Patent space carve-out | Constrained de novo | RL + structural constraints | Activity + novelty |
| Hit-to-lead | R-group replacement | TL on lead + RL | Activity + Lipinski |
| ADMET-aware design | De novo or optimization | RL | hERG + CYP + AMES + QED |

## REINVENT 4 Setup

REINVENT 4 uses a TOML configuration file specifying generator, algorithm, prior model, and scoring functions.

**Goal:** Configure a reinforcement-learning REINVENT 4 run with a prior, agent, sampling parameters, and a QED scoring component.

**Approach:** Build a release-matched REINVENT 4 staged-learning TOML config with `[parameters]` for the prior/agent checkpoints, `[[stage]]` blocks, and one or more `[[stage.scoring.component]]` blocks. Validate the config with the installed release because component parameters evolve between versions.

```toml
run_type = "staged_learning"
device = "cuda:0"

[parameters]
prior_file = "priors/reinvent.prior"
agent_file = "priors/reinvent.prior"
batch_size = 64
unique_sequences = true

[[stage]]
termination = "simple"
min_steps = 25
max_steps = 500

[stage.scoring]
type = "geometric_mean"

[[stage.scoring.component]]
[stage.scoring.component.QED]

[[stage.scoring.component.QED.endpoint]]
name = "QED"
weight = 1
```

```bash
# The REINVENT 4 CLI binary is `reinvent` (not `reinvent4`).
reinvent -l logfile.log config.toml
```

Output: a live stage CSV using `summary_csv_prefix`, plus the configured `chkpt_file` at stage termination or graceful interruption. Post-process the CSV to select molecules; REINVENT does not emit a checkpoint and SMILES file at every iteration by default.

## Scoring Function Design (Most Important Part)

A good scoring function:
- Returns 0-1 (normalized)
- Combines multiple endpoints
- Penalizes pathological generations (PAINS, unstable, unsynthesizable)

**Goal:** Build a multi-component generative reward that balances predicted activity, drug-likeness, synthesizability, and novelty.

**Approach:** Combine a QSAR sigmoid on pIC50, QED, SA-score reverse-sigmoid, and Tanimoto-similarity reverse-sigmoid via geometric mean so any zero component zeroes the total.

In REINVENT 4, define these under the active stage's `[stage.scoring]` section, with each component using the exact component and endpoint tables from the installed release's configuration examples. Do not reuse REINVENT 3 `[scoring_function]` or `[[scoring_function.components]]` syntax in a REINVENT 4 config. The accompanying example is deliberately limited to built-in, documented component structure; add predictive-property endpoints only after validating their release-specific model-container parameters.

`geometric_mean` ensures all components must be reasonably high (one zero -> zero total). `arithmetic_mean` allows compensation.

## Multi-Parameter Optimization (MPO)

Lead optimization commonly involves multiple objectives. The following component types illustrate a project-specific scoring design; weights and transforms must be fit to the actual assays and decision context.

| Component | Weight | Transformation |
|-----------|--------|----------------|
| Target activity (predicted pIC50) | 0.3 | sigmoid 5-8 |
| Selectivity (off-target ratio) | 0.2 | sigmoid 1-100 |
| QED | 0.1 | identity |
| Synthetic accessibility (SA score) | 0.1 | reverse sigmoid 1-4 |
| hERG predicted prob | 0.1 | reverse sigmoid 0.3-0.7 |
| AMES predicted prob | 0.1 | reverse sigmoid 0.3-0.7 |
| Tanimoto novelty vs known | 0.1 | reverse sigmoid 0.4-0.6 |

The weights and transformation bounds above are repository starting examples only. Normalize weights as required by the selected aggregation and tune every bound against project assay distributions and prospective behavior.

## Reward Hacking (Production Pitfall)

RL agents will find ways to maximize reward without learning the intended behavior:
- Trivial scaffolds that score high on QED
- Repeat structural motifs that game similarity scoring
- Out-of-distribution molecules that exploit QSAR overconfidence
- Trivial SMILES (e.g., "CCC...C") that match generic scoring

**Mitigations:**
- Include one or more synthesis-aware signals when they have been validated for the project; SA score alone does not establish route feasibility
- Use ensemble QSAR with uncertainty (penalize high-uncertainty predictions)
- Include diversity bonus (Tanimoto to reference)
- Add fingerprint similarity penalty within batch (prevent mode collapse)
- Validate generated samples on held-out QSAR test set

## Synthetic Accessibility Scoring

`sa_score` (Ertl 2009) measures synthetic accessibility: 1 (easy) to 10 (very hard).

```python
from rdkit.Contrib.SA_Score import sascorer
from rdkit import Chem

def sa_score(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    return sascorer.calculateScore(mol)
```

`sascorer` is shipped in RDKit Contrib in current RDKit distributions. Use the namespaced import above; do not install an unrelated top-level package.

The score is a 1-to-10 heuristic derived from fragment contributions and molecular complexity, with lower values intended to indicate easier synthesis. It is not a route planner, cost estimate, or calibrated feasibility probability. Use it as one audited reward component or annotation, never as an absolute filter.

## Diffusion-Based Generation (Modern Alternatives)

| Tool | Approach | Strength | Status |
|------|----------|----------|--------|
| DiGress (Vignac 2023) | Discrete diffusion on graphs | Conditional generation | Public |
| DiffSMol (Chen 2025) | Equivariant diffusion | 3D molecule generation | Public |
| MolDiff (Peng 2023) | Full-atom diffusion | Joint atom/bond generation | Public |
| TargetDiff (Guan 2023) | Pocket-conditioned equivariant diffusion | Structure-based design | Public |

Diffusion models iteratively denoise molecular representations, whereas REINVENT generators autoregressively construct SMILES. Diversity, validity, and drug-likeness depend on the model, training data, conditioning, and evaluation protocol; compare them on a matched benchmark for the intended task.

## Constrained / Goal-Directed Generation

**Goal:** Enforce hard structural requirements (e.g., must contain hydroxyl) and exclude PAINS without letting constraint satisfaction game the reward.

**Approach:** Stage transfer learning then RL. In REINVENT 4, `CustomAlerts` is a global structural-alert filter: a match produces zero and it is applied before score aggregation. `MatchingSubstructure` is a scoring component (1 for a match and 0.5 otherwise), so it is a soft penalty rather than a hard inclusion constraint. Apply a separate post-generation SMARTS validation step when presence of a feature is mandatory.

```toml
[[stage.scoring.component]]
[stage.scoring.component.CustomAlerts]

[[stage.scoring.component.CustomAlerts.endpoint]]
name = "Unwanted SMARTS"
params.smarts = ["PAINS_SMARTS_1", "BRENK_SMARTS_1"]

[[stage.scoring.component]]
[stage.scoring.component.MatchingSubstructure]

[[stage.scoring.component.MatchingSubstructure.endpoint]]
name = "Hydroxyl preference"
weight = 0.1
params.smarts = ["[OX2H]"]
```

There is no REINVENT 4 `filter_only` option for these components. Treat structural alerts as triage flags where appropriate, and separately verify any true hard inclusion or exclusion rule on the generated structures.

## MolMIM (NVIDIA BioNeMo)

MolMIM encodes SMILES into a learned latent space, uses gradient-free CMA-ES to optimize a user-defined property objective, and decodes candidate molecules.

```python
# Pseudo-code; requires NVIDIA NIM access
# from bionemo.molmim import MolMIMOptimizer
# optimizer = MolMIMOptimizer(model="molmim-property-optimizer")
# optimized = optimizer.optimize(seed_smiles, target_property="logp", target_value=2.0)
```

Tradeoffs against REINVENT depend on the oracle budget, objective, and implementation; benchmark both under matched constraints when selecting a generator.

## Per-Tool Failure Modes

### REINVENT RL -- mode collapse

**Trigger:** The reward, learning strategy, or diversity control favors a narrow chemotype.

**Mechanism:** Agent finds a high-scoring local maximum and stops exploring.

**Symptom:** Diversity and scaffold coverage collapse relative to a project-defined baseline while reward continues to rise.

**Fix:** Add diversity bonus to scoring; reduce sigma; reset agent if collapsed.

### REINVENT TL -- overfitting

**Trigger:** Transfer learning data are too small or homogeneous for the intended generalization task.

**Mechanism:** Generator memorizes training set; no generalization.

**Symptom:** Generated molecules near-identical to training set actives.

**Fix:** Use larger training set; mix with diverse external sample; apply RL after TL.

### Generated molecule unsynthesizable

**Trigger:** SA score missing from reward.

**Mechanism:** The reward omits synthesis evidence, allowing candidates with no plausible validated route to score well.

**Symptom:** AiZynthFinder cannot solve route; medchem rejects.

**Fix:** Combine audited synthesis-aware annotations with reaction- or route-based validation for selected candidates.

### PAINS in generation

**Trigger:** No structural alerts in scoring.

**Mechanism:** Curcumin / rhodanine / quinone scaffolds optimize for activity (false positives in training data).

**Symptom:** Generated molecules match PAINS_A.

**Fix:** Flag relevant structural alerts for orthogonal assay review or apply a project-justified penalty; never reward a PAINS match.

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
| Speed | Implementation-, hardware-, and oracle-dependent | Implementation-, hardware-, and sampling-step-dependent |
| Output diversity | Model/training/reward-dependent | Model/training/conditioning-dependent |
| Drug-likeness of output | Training- and reward-dependent | Training- and conditioning-dependent |
| Scoring flexibility | Excellent (TOML config) | Method-specific |
| Production maturity | High | Emerging |
| When to use | When its priors and scoring system validate well for the task | When a matched benchmark supports the selected diffusion model |

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| REINVENT generates invalid SMILES | Prior/tokenization/model mismatch or sampling issue | Inspect invalid-token logs, prior compatibility, and release-matched sampling settings; sigma is not a token sampling rate |
| QSAR score all 0.0 | Out-of-domain molecules | Ensemble + uncertainty; reject high-uncertainty |
| All generations duplicates | `unique_sequences=False` | Set `unique_sequences=true` |
| Generated SMILES too long | Prior/model sequence behavior | Use a release-documented sampling constraint or post-generation project rule; do not invent a staged-learning `max_length` field |
| Reward stuck at 0.5 | Constraints conflict | Inspect scoring components; reduce constraint count |
| Diffusion model crashes | Input violates model-specific pocket/size contract | Follow that model release's documented preprocessing and limits |
| MolMIM cold-start slow | Latent search exhaustiveness | Reduce search budget |
| Optimization converges trivially | Reward gradient dominated by one term | Use geometric_mean; rebalance weights |

## References

- Loeffler et al., *J. Cheminformatics* 16:20 (2024) -- REINVENT 4 framework and four generator families (DOI 10.1186/s13321-024-00812-5).
- Olivecrona M et al., *J. Cheminformatics* 9:48 (2017) -- REINVENT original (DOI 10.1186/s13321-017-0235-x).
- Vignac et al., *ICLR* (2023) -- DiGress discrete diffusion.
- Chen H et al., *Nat. Mach. Intell.* 7:758-770 (2025) -- DiffSMol structure-based 3D molecular generation (DOI 10.1038/s42256-025-01030-w).
- Peng X, Guan J, Liu Q, Ma J. *Proc. ICML*, PMLR 202:27611-27629 (2023) -- MolDiff full-atom molecular diffusion.
- Guan J et al., *ICLR* (2023) -- TargetDiff pocket-conditioned 3D equivariant diffusion (OpenReview: kJqXEPXMsE0).
- Reidenbach D, Livne M, Ilango RK, Gill M, Israeli J. *MLDD Workshop at ICLR* (2023) -- MolMIM and CMA-ES latent-space optimization (OpenReview: iOJlwUTUyrN).
- Jin W, Barzilay R, Jaakkola T. *Proc. ICML*, PMLR 80:2323-2332 (2018) -- JT-VAE junction-tree.
- Ertl P, Schuffenhauer A. *J. Cheminformatics* 1:8 (2009) -- SA score (DOI 10.1186/1758-2946-1-8).
- REINVENT 4 official repository and release-matched configs: https://github.com/MolecularAI/REINVENT4
- RDKit SA-score implementation: https://github.com/rdkit/rdkit/tree/master/Contrib/SA_Score

## Related Skills

- chemoinformatics/qsar-modeling - Build scoring models for generative
- chemoinformatics/retrosynthesis - Validate synthetic feasibility post-generation
- chemoinformatics/molecular-standardization - Standardize generated SMILES
- chemoinformatics/admet-prediction - ADMET in scoring components
- chemoinformatics/substructure-search - PAINS / BRENK filter for generation
- chemoinformatics/scaffold-analysis - Scaffold-aware generation control
- chemoinformatics/reaction-enumeration - Alternative to generative for combinatorial
- chemoinformatics/virtual-screening - Validate generated against target

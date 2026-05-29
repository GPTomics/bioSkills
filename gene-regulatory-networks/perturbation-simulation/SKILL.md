---
name: bio-gene-regulatory-networks-perturbation-simulation
description: Simulate transcription factor perturbation effects on cell state in silico with CellOracle and Dynamo, and predict transcriptional responses to genetic perturbations with GEARS, scGen, and CPA. Covers the direction-not-magnitude principle, local-linear validity, the GRN/velocity error it inherits, baseline discipline (mean and additive baselines), and the validation gap. Use when predicting TF knockout or overexpression effects, ranking driver TFs for fate transitions, or planning perturbation experiments. For GRN construction see multiomics-grn; for experimental Perturb-seq see single-cell/perturb-seq.
tool_type: python
primary_tool: CellOracle
---

## Version Compatibility

Reference examples tested with: CellOracle 0.18+, scanpy 1.10+, anndata 0.10+; Dynamo (dynamo-release 1.4+), GEARS, CPA (scvi-tools ecosystem) where used.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Perturbation Simulation

**"Predict what happens if I knock out this transcription factor"** -> Propagate a forced expression change through a learned GRN (or a learned velocity vector field), then project the shift onto the cell-state manifold to predict the DIRECTION cells move.
- Python: `celloracle.Oracle` for GRN-based in silico KO/overexpression
- Python: `dynamo` for vector-field (Jacobian) perturbation; `GEARS`/`CPA` for response prediction

## The Single Most Important Modern Insight -- These Methods Predict a Direction, Locally, and Inherit Every Error of the Underlying Model

In silico perturbation outputs a **shift vector whose orientation is the deliverable**; its absolute magnitude is uncalibrated. CellOracle propagates an expression *shift* (not absolute counts) and is explicit that it is for **hypothesis generation**, predicting direction not magnitude. The machinery is a **local linear approximation**: CellOracle's per-cluster GRN is a regularized linear model iterated a small fixed number of steps (`n_propagation` default 3), and Dynamo's perturbation is a first-order Taylor expansion of the learned vector field via its Jacobian. So predictions are valid only for **small perturbations near the observed manifold, for one or a few steps** -- they cannot capture long-range feedback, trans effects, or transitions to attractors not already in the data (CellOracle moves probability mass among observed states on the KNN graph; it cannot invent a new cell type). And every prediction **inherits the errors of its substrate**: CellOracle is only as good as its GRN, Dynamo only as good as its RNA-velocity estimate (which has documented reliability problems, Bergen 2021; Gorin 2022).

The 2025 baseline discipline is non-negotiable: for perturbation *response prediction*, deep and foundation models do **not** consistently beat trivial baselines -- Ahlmann-Eltze, Huber & Anders 2025 (*Nat Methods* 22:1657) show that for unseen single perturbations no DL model beats predicting the training **mean** response, and for combinations none beats the **additive** baseline (summing single-gene effects). Always report the mean baseline (unseen singles) and additive baseline (combinations); a method that does not beat them adds no value. And ask the validation question of any result: how many predictions were tested against a real knockout, and were the failures reported?

## Method Taxonomy

| Method | Citation | Mechanism | Predicts | Inherits errors of |
|--------|----------|-----------|----------|--------------------|
| CellOracle | Kamimoto 2023 *Nature* | shift propagated through per-cluster linear GRN, projected on KNN graph | direction of cell-state shift | the GRN (base + regression) |
| Dynamo | Qiu 2022 *Cell* | Jacobian of a learned analytical vector field (Df = J*Dx) | redirected fate / least-action path | RNA-velocity estimate |
| GEARS | Roohani 2024 *Nat Biotechnol* | GNN over gene + GO graphs | unseen single & combinatorial responses, GI type | training Perturb-seq + graphs |
| scGen / CPA | Lotfollahi 2019/2023 | latent-space arithmetic / disentangled composable embeddings | response to known/composed perturbations | training distribution |
| Boolean / ODE (BoolNet, BoolODE) | Müssel 2010 | mechanistic logic/kinetic simulation; attractors | extrapolative fate, multistability | the hand-curated wiring |

Data-driven methods (CellOracle/Dynamo/DL) are genome-scale, local, linear, direction-only; mechanistic Boolean/ODE models are hand-curated and small but nonlinear and extrapolative. State which regime a result belongs to.

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Predict TF KO/OE effect on cell fate from scRNA + accessibility prior | CellOracle | GRN-based shift onto the developmental flow; direction + perturbation score |
| Have high-quality (labeling-based) velocity, want fate redirection | Dynamo | learned vector field + Jacobian; no GRN prior needed |
| Predict response to unseen perturbations / combinations | GEARS or CPA | generalize via gene/GO graphs or composable embeddings -- but check baselines |
| Need mechanistic attractor/multistability reasoning | Boolean/ODE (BoolNet) | only family that extrapolates beyond observed states |
| Build the base GRN that CellOracle perturbs | -> multiomics-grn | base GRN construction lives there |
| Validate predictions experimentally | -> single-cell/perturb-seq | Perturb-seq is the interventional ground truth |

## CellOracle: Build the GRN

**Goal:** Learn context-specific (per-cluster) regulatory weights on top of a base-GRN prior.

**Approach:** Import scRNA-seq and the base GRN, impute, then fit a regularized linear model per cluster and filter links by significance.

```python
import scanpy as sc
import celloracle as co
import pandas as pd

adata = sc.read_h5ad('clustered.h5ad')          # normalized, log, PCA, clustered
oracle = co.Oracle()
oracle.import_anndata_as_raw_count(adata=adata, cluster_column_name='cell_type',
                                   embedding_name='X_umap')
oracle.import_TF_data(TF_info_matrix=pd.read_parquet('base_grn.parquet'))  # prior; see multiomics-grn
oracle.perform_PCA()
oracle.knn_imputation(n_pca_dims=50, k=None, balanced=True, b_sight=3000, b_maxl=1500)

links = oracle.get_links(cluster_name_for_GRN_unit='cell_type', alpha=10)
links.filter_links(p=0.001, weight='coef_abs', threshold_number=2000)
oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
oracle.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFdict=True)
```

## CellOracle: Simulate and Score

**Goal:** Predict the direction cells move under a TF knockout or overexpression.

**Approach:** Set the gene of interest to 0 (KO) or above-max (OE), propagate the shift a few steps, estimate KNN transition probabilities, compute the per-cell embedding shift, and score it against the developmental vector field by inner product.

```python
import numpy as np

# perturb_condition sets the clamped value: 0.0 = knockout; above observed max = overexpression.
# n_propagation is small BY DESIGN (3); increasing it amplifies linear-approximation error.
oracle.simulate_shift(perturb_condition={'GATA1': 0.0}, n_propagation=3)
oracle.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1)
oracle.calculate_embedding_shift(sigma_corr=0.05)

# Perturbation score = inner product of the simulated shift with the developmental flow.
# It is meaningless without a defined development vector field (Gradient_calculator below).
from celloracle.applications import Gradient_calculator
grad = Gradient_calculator(oracle_object=oracle, pseudotime_key='pseudotime')
grad.calculate_p_mass(smooth=0.8, n_grid=40, n_neighbors=200)
grad.calculate_mass_filter(min_mass=0.01, plot=False)   # required before calculate_gradient
grad.calculate_gradient()
shift = np.sqrt((oracle.adata.obsm['delta_embedding'] ** 2).sum(axis=1))   # magnitude is uncalibrated
```

## Dynamo: Vector-Field Perturbation

**Goal:** Predict fate redirection from a learned dynamical system rather than a GRN prior.

**Approach:** Reconstruct the analytical vector field, compute its Jacobian, then apply a genetic perturbation as a local linear response.

```python
import dynamo as dyn

dyn.vf.VectorField(adata, basis='umap')          # learn the analytical vector field first
dyn.vf.jacobian(adata, regulators=['GATA1'], effectors=['SPI1'])
# expression is a LIST aligned to genes (negative = suppress); perturbation writes a new basis.
dyn.pd.perturbation(adata, 'GATA1', [-100], emb_basis='umap')
dyn.pl.streamline_plot(adata, color='cell_type', basis='umap_perturbation')
```

## Baseline Discipline (always run these)

**Goal:** Establish whether a perturbation-prediction model beats triviality before trusting it.

**Approach:** Compare predictions against the training-mean response (unseen single perturbations) and the additive baseline (combinations).

```python
import numpy as np

# Mean baseline: predict the average perturbed profile across the training perturbations.
mean_baseline = train_perturbed.mean(axis=0)

# Additive baseline for a double perturbation: sum the two single-gene log-fold-changes.
additive_pred = lfc_singleA + lfc_singleB

# A model only "works" if its error is below these baselines (Ahlmann-Eltze 2025).
```

## Per-Method Failure Modes

### Reporting magnitudes as quantitative
**Trigger:** citing a predicted fold-change as an effect size. **Mechanism:** the methods are direction-only; amplitude is uncalibrated. **Symptom:** quantitative KO transcriptome claims. **Fix:** report direction / perturbation score; validate magnitude experimentally.

### No trivial baseline
**Trigger:** a DL/foundation-model perturbation predictor reported without baselines. **Mechanism:** mean (singles) and additive (combos) baselines are often competitive. **Symptom:** "state of the art" with no mean/additive comparison. **Fix:** report both baselines; require the model to beat them.

### Perturbing far outside the manifold
**Trigger:** extreme overexpression, knocking out a gene not expressed in the cluster, or expecting a brand-new cell type. **Mechanism:** local linear validity is violated; CellOracle can only move mass among observed states. **Symptom:** implausible jumps. **Fix:** keep perturbations small and near observed states; interpret one step.

### Cranking n_propagation
**Trigger:** raising `n_propagation` to get "long-range" effects. **Mechanism:** iterating a linear model amplifies approximation error. **Symptom:** runaway shifts. **Fix:** keep the small default; treat it as a local estimate.

### Trusting the substrate blindly
**Trigger:** treating the CellOracle GRN as ground truth, or running Dynamo on noisy splicing-based velocity. **Mechanism:** predictions inherit GRN/velocity error (Bergen 2021; Gorin 2022; and a 2024 preprint critiques CellOracle's GRN for ignoring distal interactions). **Symptom:** confident predictions on a poorly-validated network/field. **Fix:** sanity-check the GRN/velocity; present results as hypotheses.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| n_propagation = 3 | CellOracle default | small by design; local linear validity |
| KO value = 0.0; OE = above observed max | CellOracle convention | clamped perturbation input |
| link filter p<0.001, top ~2000 by coef_abs | CellOracle tutorial | sparsify the GRN before simulation |
| mean baseline (unseen singles), additive (combos) | Ahlmann-Eltze 2025 | the bar any predictor must clear |
| sigma_corr ~0.05 | CellOracle default | embedding-shift smoothing |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| perturbation score all ~0 / meaningless | no development vector field defined | run `Gradient_calculator` before scoring |
| KeyError on gene in `simulate_shift` | gene absent from the fitted GRN | confirm the gene is in the base GRN and expressed |
| Dynamo perturbation errors | vector field/Jacobian not computed | run `dyn.vf.VectorField` then `dyn.vf.jacobian` first |
| huge predicted shift | n_propagation too high / out-of-manifold perturbation | keep n_propagation small; perturb near observed states |
| "model beats prior methods" but no baseline | missing mean/additive comparison | add baselines (Ahlmann-Eltze 2025) |

## References

- Kamimoto K, et al. 2023. Dissecting cell identity via network inference and in silico gene perturbation (CellOracle). *Nature* 614(7949):742-751.
- Qiu X, et al. 2022. Mapping transcriptomic vector fields of single cells (Dynamo). *Cell* 185(4):690-711.e45.
- Roohani Y, Huang K, Leskovec J. 2024. Predicting transcriptional outcomes of novel multigene perturbations with GEARS. *Nat Biotechnol* 42(6):927-935.
- Lotfollahi M, Wolf FA, Theis FJ. 2019. scGen predicts single-cell perturbation responses. *Nat Methods* 16:715-721.
- Lotfollahi M, et al. 2023. Predicting cellular responses to complex perturbations in high-throughput screens (CPA). *Mol Syst Biol* 19(6):e11517.
- Ahlmann-Eltze C, Huber W, Anders S. 2025. Deep-learning-based gene perturbation effect prediction does not yet outperform simple linear baselines. *Nat Methods* 22:1657-1661.
- Bergen V, et al. 2021. RNA velocity - current challenges and future perspectives. *Mol Syst Biol* 17(8):e10282.
- Gorin G, Fang M, Chari T, Pachter L. 2022. RNA velocity unraveled. *PLoS Comput Biol* 18(9):e1010492.
- Müssel C, Hopfensitz M, Kestler HA. 2010. BoolNet. *Bioinformatics* 26(10):1378-1380.

## Related Skills

- multiomics-grn - build the CellOracle base GRN from accessibility
- scenic-regulons - regulon activity as an alternative TF-driver readout
- single-cell/perturb-seq - experimental Perturb-seq, the interventional ground truth
- single-cell/trajectory-inference - pseudotime/development flow for the perturbation score
- causal-genomics/mendelian-randomization - population-genetics route to causality (contrast)

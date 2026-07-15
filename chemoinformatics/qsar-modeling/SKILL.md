---
name: bio-qsar-modeling
description: Builds QSAR / QSPR models using chemprop D-MPNN, MolFormer, Uni-Mol, ChemBERTa, random forest baselines, and Gaussian processes with explicit handling of OECD 5 principles, applicability domain (kNN, leverage, conformal prediction, Mahalanobis), scaffold-balanced splits, ensemble uncertainty, calibration (Platt, isotonic), feature importance (SHAP, atomic attribution), and prospective validation. Use when building target-specific predictive models from in-house bioassay data, ADMET endpoints, or selectivity profiles.
tool_type: python
primary_tool: chemprop
---

## Version Compatibility

Reference examples target: chemprop 2.2.x (major API change from 1.x), RDKit 2024.09+, scikit-learn >=1.4,<1.6, MAPIE >=0.8,<1.0 for the `MapieRegressor` example, shap 0.44+, and pytorch 2.1+. Recheck examples before widening these bounds because Chemprop, scikit-learn calibration, and MAPIE interfaces evolve independently.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `chemprop train --help` (chemprop 2.x); `chemprop_train --help` (1.x legacy)

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# QSAR Modeling

Build quantitative structure-activity relationship models from molecular structure inputs. The choice of model, featurization, and split strategy determines whether the model captures transferable chemical signal or memorizes the training data. chemprop D-MPNN with optional Morgan / RDKit descriptors is a useful open-source approach; transformer-based methods (MolFormer, Uni-Mol, ChemBERTa) should be compared on the same split and endpoint. The OECD validation principles support transparent documentation and evaluation of (Q)SAR models, but following them does not by itself confer regulatory acceptance.

For descriptor/fingerprint choices, see `chemoinformatics/molecular-descriptors`. For ADMET-specific QSAR, see `chemoinformatics/admet-prediction`. For molecular standardization (critical upstream), see `chemoinformatics/molecular-standardization`.

## Model Taxonomy

| Model | Architecture | Use case | Fails when |
|-------|--------------|----------|------------|
| Random Forest + ECFP4 | Classical baseline | Small-data comparison, interpretability | May miss signal not represented by the fingerprint |
| chemprop D-MPNN | Directed message passing | Graph-learning candidate to benchmark | Can overfit when data are sparse or biased |
| chemprop D-MPNN + RDKit 2D | Hybrid graph + descriptors | Useful hybrid baseline; compare on the same split | Diminishing returns at large data |
| MolFormer | SMILES transformer | Large public training data benefit | Compute overhead; OOD risk |
| Uni-Mol | 3D-aware transformer | 3D-relevant endpoints (binding) | Requires 3D conformers |
| ChemBERTa-2 | SMILES transformer pretrained on up to 77M molecules | SMILES language-model baseline | Fine-tuning benefit is endpoint- and split-dependent |
| Gaussian Process + ECFP4 | Probabilistic | Active learning; uncertainty | O(N^3) scaling |
| MultiTask DNN | Joint training | Multiple endpoints | Data must overlap |

**Decision:** Compare a fingerprint-based baseline with chemprop under the same split and endpoint. Add a pretrained transformer or 3D model only when its representation, compute cost, and validation design fit the deployment question; dataset size alone does not determine the winner.

## Decision Tree by Scenario

| Dataset context | Endpoint type | Model to benchmark |
|--------------|---------------|-------|
| Sparse labels or few independent series | Regression / classification | Regularized fingerprint baseline; quantify instability and avoid unsupported deployment |
| Multiple scaffold groups with adequate labels | Regression / classification | Fingerprint baseline plus chemprop on identical splits |
| Large public or internal training collection | Regression / classification | Benchmark chemprop and a relevant pretrained representation |
| Multi-task | Related endpoints (CYP3A4, CYP2D6, etc.) | chemprop MultiTask |
| 3D-relevant | Binding, conformer-dependent | Uni-Mol with conformer ensemble |

## OECD 5 Principles

The OECD principles were agreed in 2004; the 2007 guidance explains their application:

1. **Defined endpoint**: specific bioassay, units, threshold definitions
2. **Unambiguous algorithm**: reproducible code, fixed random seeds, version-pinned dependencies
3. **Defined applicability domain (AD)**: where the model is valid
4. **Appropriate measures of goodness-of-fit, robustness, and predictivity**: external test set and suitable validation
5. **Mechanistic interpretation, if possible**: biological/chemical rationale where available

For non-regulatory QSAR, all 5 still good practice; especially **AD definition** is critical.

## Applicability Domain Methods

| Method | Definition | Pro | Con |
|--------|-----------|-----|-----|
| **Ensemble variance** | Std across N-model ensemble predictions | Supported by `chemprop predict --uncertainty-method ensemble` when multiple model paths are supplied | Assumes useful ensemble diversity; not calibrated coverage |
| kNN distance | Mean Tanimoto to k nearest in training | Easy to interpret | Doesn't account for label distribution |
| Leverage | Hat matrix diagonal | Statistical | Linear assumptions |
| KDE on PCA | Density in feature space | Captures multivariate structure | Density choice subjective |
| Mahalanobis distance | Covariance-aware distance | Theoretically motivated | High-dim instability |
| Conformal prediction | Per-prediction interval or set | Finite-sample marginal coverage under exchangeability | Requires a calibration design and compatible predictor |
| Bayesian / MC-dropout | Posterior or dropout variance | Direct uncertainty | Computational cost |
| Tanimoto coverage | At least 1 NN within threshold | Practical | Threshold subjective |

Ensemble disagreement is one useful uncertainty diagnostic, not a formally defined applicability domain or calibrated coverage guarantee. If using a threshold such as a training-distribution percentile, label it as a project-defined heuristic and validate it prospectively.

## chemprop 2.x Training (CLI)

**Goal:** Train five replicated chemprop runs, each containing a five-model D-MPNN ensemble with RDKit 2D descriptor features and a scaffold-balanced train/validation/test split.

**Approach:** For current chemprop 2.x, invoke `chemprop train` with `--molecule-featurizers rdkit_2d`, `--num-replicates 5`, `--ensemble-size 5`, and `--split scaffold_balanced`. Replicates repeat splitting/training with incremented seeds; they are not five-fold cross-validation. Confirm the exact flags with `chemprop train --help` because the v2 CLI continues to evolve.

```bash
# chemprop 2.x CLI (current): use 'chemprop train' (space; dashes not underscores)
chemprop train \
    --data-path data.csv \
    --task-type classification \
    --save-dir model_dir \
    --molecule-featurizers rdkit_2d \
    --num-replicates 5 \
    --ensemble-size 5 \
    --epochs 50 \
    --batch-size 128 \
    --split scaffold_balanced \
    --split-sizes 0.8 0.1 0.1 \
    --metric roc

# chemprop 1.x legacy CLI (for backwards reference):
# chemprop_train --data_path data.csv --dataset_type classification ...
```

Key flags (chemprop 2.x):
- `--molecule-featurizers rdkit_2d`: include current v2 RDKit descriptors, which are scaled by default (the legacy v1-normalized generator is `v1_rdkit_2d_normalized`)
- `--num-replicates 5`: repeat the split/training workflow with successive seeds; this replaced `--num-folds` in chemprop 2.1
- `--ensemble-size 5`: train five models per replicate for an ensemble prediction
- `--split scaffold_balanced`: prevent scaffold leakage (was `--split_type` in 1.x)
- `--split-sizes 0.8 0.1 0.1`: 80/10/10 train/val/test

Total models: 25 (5 replicates x 5 ensemble members). Report which predictions are being aggregated and treat ensemble standard deviation as an uncertainty diagnostic, not a calibrated guarantee.

At prediction time, uncertainty output is opt-in and requires the actual saved model paths:

```bash
chemprop predict --test-path test.csv \
    --model-paths path/to/model_1.ckpt path/to/model_2.ckpt \
    --uncertainty-method ensemble \
    --preds-path predictions.csv
```

## Scaffold-Balanced Split

**Goal:** Partition a SMILES dataset into train/val/test such that no Bemis-Murcko scaffold appears in more than one split (prevents chemotype leakage).

**Approach:** Group compounds by scaffold and assign whole scaffold groups to train, validation, or test. chemprop's `--split scaffold_balanced` implements a scaffold-based allocation; `--class-balance` is a separate training option and does not make this split outcome-stratified. The chemprop default split is random, so request scaffold-balanced explicitly when it matches the deployment question.

`scaffold_balanced` assigns each scaffold group to one of train / validation / test, reducing direct scaffold leakage. It is not universally the correct validation design: time splits, externally defined series, grouped cross-validation, and prospective tests may better represent a particular deployment setting.

Choose and document the primary split before model selection. A random split can answer an interpolation question but often shares close analogues across partitions; a scaffold split tests transfer across scaffold groups; a time or prospective split tests the historical deployment process. If several splits are reported, interpret their differences as split-specific sensitivity rather than a universal "true generalization gap."

## Conformal Prediction for Calibrated Uncertainty

Use conformal prediction when calibrated marginal coverage under the stated exchangeability assumptions matters. Ensemble variance is simpler, but it is not a substitute for a conformal guarantee.

```python
# MAPIE expects a scikit-learn-compatible estimator (.fit / .predict / .predict_proba).
# chemprop 2.x is NOT scikit-learn-compatible out of the box -- either wrap chemprop
# in a thin sklearn estimator class or use MAPIE only with the sklearn baseline.
from mapie.regression import MapieRegressor
from sklearn.ensemble import RandomForestRegressor

base = RandomForestRegressor(n_estimators=500, random_state=42)
mapie = MapieRegressor(estimator=base, method='plus', cv=5)
mapie.fit(X_train, y_train)
y_pred, y_intervals = mapie.predict(X_test, alpha=0.1)  # alpha=0.1 -> 90% coverage
```

Alpha 0.05 targets 95% marginal coverage and alpha 0.10 targets 90%, subject to the conformal method's assumptions. MAPIE supports the sklearn baseline directly; integrating chemprop requires a separately implemented and tested compatible wrapper.

## SHAP / Atomic Attribution

For mechanistic interpretation:

For a scikit-learn-style model (e.g., Random Forest baseline on ECFP4), SHAP integrates directly:

```python
import shap
from sklearn.ensemble import RandomForestClassifier

# X_train / X_test are Morgan fingerprint arrays (n_samples, n_bits)
model = RandomForestClassifier(n_estimators=500, random_state=42).fit(X_train, y_train)
explainer = shap.TreeExplainer(model)
shap_values = explainer.shap_values(X_test)

# Per-bit contribution; for atomic interpretation, map bits back to
# generating atoms via AllChem.GetMorganFingerprintAsBitVect(mol, ..., bitInfo=bi)
# and aggregate SHAP across all bits triggered by each atom.
```

For chemprop D-MPNN, SHAP requires a custom wrapper (chemprop is not sklearn-compatible). A PyTorch attribution method must be adapted to the model's graph inputs and validated; the chemprop 2.x CLI does not provide the `--uncertainty-method classification` atom-attribution interface. Use directly supported fingerprint SHAP for the classical baseline unless a tested graph-attribution implementation is available.

## Bayesian Optimization for Active Learning

```python
import numpy as np
from scipy.stats import norm
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF

gp = GaussianProcessRegressor(kernel=RBF(length_scale=1.0), random_state=42)
gp.fit(X_train, y_train)
mu, sigma = gp.predict(X_pool, return_std=True)

# Expected Improvement
def expected_improvement(mu, sigma, y_best, xi=0.01):
    improvement = mu - y_best - xi
    ei = np.zeros_like(mu, dtype=float)
    nonzero = sigma > 0
    z = improvement[nonzero] / sigma[nonzero]
    ei[nonzero] = (
        improvement[nonzero] * norm.cdf(z)
        + sigma[nonzero] * norm.pdf(z)
    )
    return ei

ei = expected_improvement(mu, sigma, y_train.max())
next_to_test = X_pool[ei.argmax()]
```

For chemprop + active learning, replace GP with chemprop ensemble + ensemble variance.

## Calibration (Platt / Isotonic)

Deep learning probabilities are not guaranteed to be calibrated. Use Platt (logistic) or isotonic calibration for binary probabilities, choosing the method with a held-out calibration set. `--metric roc` evaluates ranking and does not automatically calibrate chemprop probabilities; export validation probabilities and fit the calibrator externally:

```python
from sklearn.isotonic import IsotonicRegression
iso = IsotonicRegression(out_of_bounds='clip').fit(val_chemprop_probs, val_true)
test_calibrated = iso.predict(test_chemprop_probs)
```

## Multi-Task QSAR

Train multiple related endpoints jointly:

```python
df = pd.DataFrame({
    'smiles': [...],
    'CYP1A2_inhibition': [...],
    'CYP2D6_inhibition': [...],
    'CYP3A4_inhibition': [...],
})
df.to_csv('multitask.csv', index=False)
```

```bash
chemprop train --data-path multitask.csv --task-type classification \
               --target-columns CYP1A2_inhibition CYP2D6_inhibition CYP3A4_inhibition \
               --save-dir multitask_model
```

Multitask learning can help when endpoints share predictive signal or data, but negative transfer is also possible. Compare single-task and multitask models under identical splits rather than assuming improvement from endpoint relatedness.

## Per-Tool Failure Modes

### Random split for QSAR

**Trigger:** Default sklearn `train_test_split`.

**Mechanism:** Compounds from same scaffold scatter across train/test; performance optimistic.

**Symptom:** Performance drops substantially from random splits to scaffold, time, external-series, or prospective evaluation.

**Fix:** Use `--split scaffold_balanced` in chemprop 2.x (or `--split_type scaffold_balanced` in chemprop 1.x legacy); or `scaffold_split` from `chemoinformatics/scaffold-analysis`.

### Class imbalance not handled

**Trigger:** 10:1 negative:positive ratio in dataset.

**Mechanism:** Default loss treats classes equally; model learns majority class.

**Symptom:** High accuracy but precision/recall on minority class poor.

**Fix:** Class-weighted loss; SMOTE; or report AUC/F1 not accuracy.

### Over-engineered features

**Trigger:** Including hundreds of descriptors (e.g., `rdkit_2d` not normalized).

**Mechanism:** Some descriptors dominate scaling; model overfits.

**Symptom:** Validation performance differs widely across runs; high feature importance noise.

**Fix:** In current chemprop 2.x use `rdkit_2d`, which is scaled by default, or supply a documented descriptor set with preprocessing fit only on the training data.

### Missing AD assessment

**Trigger:** Predicting on novel chemotypes without AD check.

**Mechanism:** Model extrapolates; predictions unreliable.

**Symptom:** Confident predictions but actual values different.

**Fix:** Predefine and validate one or more domain/uncertainty diagnostics, such as neighborhood similarity, ensemble disagreement, or conformal output, and report what each diagnostic does and does not guarantee.

### chemprop 1.x vs 2.x confusion

**Trigger:** Code/tutorial from before late 2024.

**Mechanism:** Major API change: `chemprop_train` -> `chemprop train`; Python API redesigned.

**Symptom:** ImportError or different keyword arguments.

**Fix:** Use `chemprop --version`; check 2.x documentation; migrate APIs.

### Pretrained Transformer overhead without data benefit

**Trigger:** Adding a pretrained transformer without a matched baseline and deployment-relevant validation.

**Mechanism:** The pretrained representation, fine-tuning design, and endpoint may not provide additional transferable signal.

**Symptom:** No improvement over chemprop; slower training.

**Fix:** Compare against fingerprint and chemprop baselines on the same split, and retain the transformer only when the measured benefit justifies its cost.

### Validation leakage via standardization

**Trigger:** Standardization rules or learned preprocessing parameters are chosen or fit using validation/test data.

**Mechanism:** Test-set information influences representations, feature selection, scaling, or deduplication decisions.

**Symptom:** Re-fitting preprocessing on training data alone reduces held-out performance or changes membership across splits.

**Fix:** Freeze chemistry rules before evaluation and fit learned preprocessing on training data only. Apply the frozen pipeline to validation, test, and prospective compounds while preserving endpoint-relevant stereochemistry.

## Reconciliation: Classical RF vs chemprop vs Transformer

| Aspect | RF + ECFP4 | chemprop D-MPNN | MolFormer |
|--------|-----------|-----------------|-----------|
| Data regime | Useful baseline across sizes; especially important in small data | Compare when graph learning is plausible | Compare when pretrained representations and compute are justified |
| Interpretability | Fingerprint importance or SHAP, with bit-to-atom mapping caveats | Graph attribution requires a custom, validated implementation | Model-specific attribution requires validation |
| Uncertainty | Bootstrap or conformal wrapper | Ensemble disagreement; calibrate separately when needed | Method-dependent; validate empirically |
| Hardware | CPU | CPU or GPU depending on scale | Usually GPU for fine-tuning |
| OOD performance | Benchmark on the intended split/domain | Benchmark on the intended split/domain | Benchmark on the intended split/domain |
| Production deployment | Version-pinned sklearn artifact or service | Version-pinned native checkpoint/service; do not assume ONNX support | Version-pinned framework artifact/service |

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| chemprop hangs at start | GPU OOM | Reduce batch_size; check CUDA |
| All predictions same value | Constant target | Standardize labels |
| AUC mismatched across folds | Random seed not set | `--seed 42` |
| Test AUC = train AUC | No held-out data | Use scaffold_balanced split |
| Ensemble variance always small | Ensemble members insufficiently diverse | Check the documented seed behavior and training randomness for each replicate/member |
| SHAP fails on D-MPNN | Graph inputs are not compatible with the tree-model interface | Use a tested graph-attribution implementation or report the fingerprint baseline attribution |
| MolFormer fine-tune slow | All parameters trained | Use LoRA or freeze early layers |
| Calibration degrades held-out results | Calibrator overfit or distribution shifted | Refit on a proper calibration split and report uncalibrated and calibrated metrics |

## References

- Yang K et al. "Analyzing Learned Molecular Representations for Property Prediction." *J. Chem. Inf. Model.* 59:3370–3388 (2019). DOI: 10.1021/acs.jcim.9b00237.
- Heid E et al. "Chemprop: A Machine Learning Package for Chemical Property Prediction." *J. Chem. Inf. Model.* 64:9–17 (2024). DOI: 10.1021/acs.jcim.3c01250.
- Wu Z et al. "MoleculeNet: a benchmark for molecular machine learning." *Chem. Sci.* 9:513–530 (2018). DOI: 10.1039/C7SC02664A.
- Ross J, Belgodere B, Chenthamarakshan V, Padhi I, Mroueh Y, Das P. "Large-scale chemical language representations capture molecular structure and properties." *Nat. Mach. Intell.* 4:1256–1264 (2022). DOI: 10.1038/s42256-022-00580-7.
- Zhou G, Gao Z, Ding Q et al. "Uni-Mol: A Universal 3D Molecular Representation Learning Framework." *ICLR* (2023). OpenReview: https://openreview.net/forum?id=6K2RM6wVqKu.
- Ahmad W, Simon E, Chithrananda S, Grand G, Ramsundar B. "ChemBERTa-2: Towards Chemical Foundation Models." arXiv:2209.01712 (2022). DOI: 10.48550/arXiv.2209.01712.
- OECD. "The OECD Principles for the Validation, for Regulatory Purposes, of (Q)SAR Models" (agreed 2004); *Guidance Document on the Validation of (Quantitative) Structure-Activity Relationship [(Q)SAR] Models*, No. 69 (2007). DOI: 10.1787/9789264085442-en.
- Cortés-Ciriano I, Bender A. "Concepts and Applications of Conformal Prediction in Computational Drug Discovery." arXiv:1908.03569 (2019). DOI: 10.48550/arXiv.1908.03569.
- Svensson F et al. "Conformal Regression for Quantitative Structure–Activity Relationship Modeling—Quantifying Prediction Uncertainty." *J. Chem. Inf. Model.* 58:1132–1140 (2018). DOI: 10.1021/acs.jcim.8b00054.
- Chemprop 2.x CLI documentation, training and prediction: https://chemprop.readthedocs.io/en/latest/tutorial/cli/.
- MAPIE 0.8 documentation for the version-bounded `MapieRegressor` interface: https://mapie.readthedocs.io/en/v0.8.6/.

## Related Skills

- chemoinformatics/molecular-descriptors - Featurization choices
- chemoinformatics/molecular-standardization - Mandatory upstream
- chemoinformatics/scaffold-analysis - Bemis-Murcko split implementation
- chemoinformatics/admet-prediction - ADMET-specific QSAR
- chemoinformatics/generative-design - QSAR as scoring component
- machine-learning/model-validation - General ML validation principles
- machine-learning/biomarker-discovery - Adjacent ML approaches

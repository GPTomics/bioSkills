---
name: bio-qsar-modeling
description: Builds QSAR / QSPR models using chemprop D-MPNN, MolFormer, Uni-Mol, ChemBERTa, random forest baselines, and Gaussian processes with explicit handling of OECD 5 principles, applicability domain (kNN, leverage, conformal prediction, Mahalanobis), scaffold-balanced splits, ensemble uncertainty, calibration (Platt, isotonic), feature importance (SHAP, atomic attribution), and prospective validation. Use when building target-specific predictive models from in-house bioassay data, ADMET endpoints, or selectivity profiles.
tool_type: python
primary_tool: chemprop
---

## Version Compatibility

Reference examples tested with: chemprop 2.0+ (major API change from 1.x), RDKit 2024.09+, scikit-learn 1.4+, MAPIE 0.8+ (conformal prediction), shap 0.44+, pytorch 2.1+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `chemprop train --help` (chemprop 2.x); `chemprop_train --help` (1.x legacy)

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# QSAR Modeling

Build quantitative structure-activity relationship models from molecular structure inputs. The choice of model + featurization + split strategy determines whether the model captures real chemical signal or memorizes the training data. chemprop D-MPNN with optional Morgan / RDKit descriptors is the modern open-source standard; transformer-based methods (MolFormer, Uni-Mol, ChemBERTa) compete on benchmarks but offer minimal practical gain when domain-specific data is sparse. The OECD 5 principles structure the model for regulatory acceptance: defined endpoint, unambiguous algorithm, defined applicability domain, validation, and mechanistic interpretation.

For descriptor/fingerprint choices, see `chemoinformatics/molecular-descriptors`. For ADMET-specific QSAR, see `chemoinformatics/admet-prediction`. For molecular standardization (critical upstream), see `chemoinformatics/molecular-standardization`.

## Model Taxonomy

| Model | Architecture | Use case | Fails when |
|-------|--------------|----------|------------|
| Random Forest + ECFP4 | Classical baseline | Small data (<200 compounds), interpretability | Saturates at ~AUC 0.85 for complex endpoints |
| chemprop D-MPNN | Directed message passing | Modern default; 100-10k compounds | Very small datasets (<100) overfits |
| chemprop D-MPNN + RDKit 2D | Hybrid graph + descriptors | hERG SOTA (Cai 2024 AUC 0.956) | Diminishing returns at large data |
| MolFormer | Transformer (87M params) | Large public training data benefit | Compute overhead; OOD risk |
| Uni-Mol | 3D-aware transformer | 3D-relevant endpoints (binding) | Requires 3D conformers |
| ChemBERTa-2 | Transformer (77M params) | SMILES language model | Large training data needed |
| Gaussian Process + ECFP4 | Probabilistic | Active learning; uncertainty | O(N^3) scaling |
| MultiTask DNN | Joint training | Multiple endpoints | Data must overlap |
| AttentiveFP | GNN with attention | Pre-chemprop SOTA | Less actively maintained |
| KPGT / Graphormer / GROVER | Pretrained GNN | Pretraining lift | Diminishing returns |

**Decision:** For 200-10k compounds, **chemprop 2.0 D-MPNN + RDKit 2D descriptors** is the modern standard. For <200 compounds, **Random Forest + ECFP4** is competitive and more interpretable.

## Decision Tree by Scenario

| Dataset size | Endpoint type | Model |
|--------------|---------------|-------|
| <50 | Anything | Don't model; use as test set for literature models |
| 50-200 | Regression / classification | RF + ECFP4 + cross-validation |
| 200-1k | Regression | chemprop D-MPNN + RDKit 2D, 5-fold |
| 1k-10k | Anything | chemprop D-MPNN + RDKit 2D ensemble |
| 10k-100k | Classification | chemprop ensemble OR MolFormer fine-tuning |
| >100k | Anything | MolFormer / Uni-Mol pretrained + LoRA fine-tune |
| Multi-task | Related endpoints (CYP3A4, CYP2D6, etc.) | chemprop MultiTask |
| 3D-relevant | Binding, conformer-dependent | Uni-Mol with conformer ensemble |

## OECD 5 Principles

For regulatory use (REACH, ECHA, FDA submissions):

1. **Defined endpoint**: specific bioassay, units, threshold definitions
2. **Unambiguous algorithm**: reproducible code, fixed random seeds, version-pinned dependencies
3. **Defined applicability domain (AD)**: where the model is valid
4. **Appropriate statistical validation**: external test set, cross-validation
5. **Mechanistic interpretation**: biological/chemical rationale where possible

For non-regulatory QSAR, all 5 still good practice; especially **AD definition** is critical.

## Applicability Domain Methods

| Method | Definition | Pro | Con |
|--------|-----------|-----|-----|
| **Ensemble variance (RECOMMENDED for chemprop)** | Std across N-model ensemble predictions | Built-in to `chemprop predict`; no extra code | Assumes ensemble diversity (random seeds give independent fits) |
| kNN distance | Mean Tanimoto to k nearest in training | Easy to interpret | Doesn't account for label distribution |
| Leverage | Hat matrix diagonal | Statistical | Linear assumptions |
| KDE on PCA | Density in feature space | Captures multivariate structure | Density choice subjective |
| Mahalanobis distance | Covariance-aware distance | Theoretically motivated | High-dim instability |
| Conformal prediction | Per-prediction confidence interval | Distribution-free guarantee | Requires calibration set + sklearn-compatible estimator |
| Bayesian / MC-dropout | Posterior or dropout variance | Direct uncertainty | Computational cost |
| Tanimoto coverage | At least 1 NN within threshold | Practical | Threshold subjective |

**Operational rule:** For chemprop, set `--ensemble-size 5` at training; at predict time, flag predictions with ensemble std > P95 of the training-set ensemble std distribution as out-of-AD. This is the standard chemprop-native AD measure and avoids the overhead of wrapping with MAPIE/conformal.

## chemprop 2.0 Training (CLI)

**Goal:** Train a 5-fold x 5-model chemprop D-MPNN ensemble with RDKit 2D descriptor features and scaffold-balanced split.

**Approach:** Invoke `chemprop train` with `--molecule-featurizers rdkit_2d_normalized`, `--num-folds 5`, `--ensemble-size 5`, and `--split scaffold_balanced` to produce 25 models for ensemble mean prediction + variance-based applicability domain.

```bash
# chemprop 2.x CLI (current): use 'chemprop train' (space; dashes not underscores)
chemprop train \
    --data-path data.csv \
    --task-type classification \
    --save-dir model_dir \
    --molecule-featurizers rdkit_2d_normalized \
    --num-folds 5 \
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
- `--molecule-featurizers rdkit_2d_normalized`: include physchem descriptors (was `--features_generator` in 1.x)
- `--num-folds 5`: 5-fold cross-validation
- `--ensemble-size 5`: 5 models per fold for uncertainty
- `--split scaffold_balanced`: prevent scaffold leakage (was `--split_type` in 1.x)
- `--split-sizes 0.8 0.1 0.1`: 80/10/10 train/val/test

Total models: 25 (5 folds x 5 ensemble). Use ensemble mean as prediction; ensemble std as uncertainty.

## Scaffold-Balanced Split (chemprop 2.0 default)

**Goal:** Partition a SMILES dataset into train/val/test such that no Bemis-Murcko scaffold appears in more than one split (prevents chemotype leakage).

**Approach:** Group compounds by scaffold and greedily assign whole scaffolds to train until target fraction, then val, then test. chemprop's `--split scaffold_balanced` implements this with optional class-stratification.

scaffold_balanced split assigns entire scaffolds to one of train / val / test, preventing chemotype leakage. Modern QSAR uses scaffold split exclusively.

| Split type | Random | Scaffold balanced |
|------------|--------|--------------------|
| Train AUC | 0.99 | 0.95 |
| Test AUC | 0.95 | 0.75-0.85 |
| Generalization | Optimistic | Realistic |

The gap between random and scaffold split is the **true generalization gap**; report both.

## Conformal Prediction for Calibrated Uncertainty

Use conformal prediction when calibrated coverage guarantees matter (regulatory submission, decision thresholds with cost asymmetry). For routine QSAR, chemprop ensemble variance is simpler and sufficient.

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

Alpha choice: 0.05 (95% coverage) for safety endpoints; 0.10 (90%) for activity ranking. For chemprop, the simpler ensemble-variance route (set `--ensemble-size 5` and use prediction std as AD signal) is preferred unless conformal coverage guarantees are specifically required.

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

For chemprop D-MPNN, SHAP requires a custom wrapper (chemprop is not sklearn-compatible). Use `shap.GradientExplainer` on the underlying PyTorch module after exposing it via a wrapper class. In practice, chemprop's built-in atom attribution (`chemprop predict --uncertainty-method classification` with atom-level outputs) is the supported route for atomic interpretation. Use SHAP only when comparing chemprop to a classical baseline.

## Bayesian Optimization for Active Learning

```python
from sklearn.gaussian_process import GaussianProcessRegressor

gp = GaussianProcessRegressor(kernel='rbf')
gp.fit(X_train, y_train)
mu, sigma = gp.predict(X_pool, return_std=True)

# Expected Improvement
def expected_improvement(mu, sigma, y_best, xi=0.01):
    z = (mu - y_best - xi) / sigma
    return (mu - y_best - xi) * norm.cdf(z) + sigma * norm.pdf(z)

ei = expected_improvement(mu, sigma, y_train.max())
next_to_test = X_pool[ei.argmax()]
```

For chemprop + active learning, replace GP with chemprop ensemble + ensemble variance.

## Calibration (Platt / Isotonic)

Deep learning native probabilities rarely well-calibrated. Apply Platt scaling on validation set:

```python
from sklearn.calibration import CalibratedClassifierCV
from sklearn.ensemble import RandomForestClassifier

# For sklearn-compatible models, calibration is direct:
base = RandomForestClassifier(n_estimators=500, random_state=42).fit(X_train, y_train)
calibrated = CalibratedClassifierCV(base, method='isotonic', cv='prefit')
calibrated.fit(X_val, y_val)
y_proba_calibrated = calibrated.predict_proba(X_test)
```

Use Platt (logistic) for binary; isotonic for general. For chemprop, native calibration is applied automatically when training with classification + `--metric roc`; for further re-calibration, save validation probabilities and fit isotonic regression externally:

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

Joint training improves performance when endpoints are correlated (Wu 2018 MoleculeNet). For independent endpoints, multitask hurts.

## Per-Tool Failure Modes

### Random split for QSAR

**Trigger:** Default sklearn `train_test_split`.

**Mechanism:** Compounds from same scaffold scatter across train/test; performance optimistic.

**Symptom:** Cross-validation AUC > 0.95; real-world performance much lower.

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

**Fix:** Use `rdkit_2d_normalized`; or restrict to <50 standard descriptors.

### Missing AD assessment

**Trigger:** Predicting on novel chemotypes without AD check.

**Mechanism:** Model extrapolates; predictions unreliable.

**Symptom:** Confident predictions but actual values different.

**Fix:** Always compute kNN Tanimoto distance + ensemble variance; flag low-confidence predictions.

### chemprop 1.x vs 2.x confusion

**Trigger:** Code/tutorial from before late 2024.

**Mechanism:** Major API change: `chemprop_train` -> `chemprop train`; Python API redesigned.

**Symptom:** ImportError or different keyword arguments.

**Fix:** Use `chemprop --version`; check 2.x documentation; migrate APIs.

### Pretrained Transformer overhead without data benefit

**Trigger:** Using MolFormer on <500 compounds.

**Mechanism:** Pretrained representation needs sufficient fine-tuning data to specialize.

**Symptom:** No improvement over chemprop; slower training.

**Fix:** For <500 compounds, stick with chemprop or RF. Use pretrained Transformers for >10k.

### Validation leakage via standardization

**Trigger:** Different standardization for train vs test.

**Mechanism:** Train compounds with R-isomer; test compounds with S-isomer not in train.

**Symptom:** Test performance better than realistic.

**Fix:** Standardize unified preprocessing; same canonicalization train + test.

## Reconciliation: Classical RF vs chemprop vs Transformer

| Aspect | RF + ECFP4 | chemprop D-MPNN | MolFormer |
|--------|-----------|-----------------|-----------|
| Data size sweet spot | 50-1k | 200-10k | >10k |
| Interpretability | High (Gini importance) | Medium (SHAP) | Low |
| Uncertainty | Bootstrap | Ensemble | MC-dropout |
| Hardware | CPU | GPU recommended | GPU mandatory |
| OOD performance | Most degraded | Better with ensemble | Worst (long-range memorization) |
| Production deployment | Easy (pickle) | OK (ONNX export) | Heavy |

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| chemprop hangs at start | GPU OOM | Reduce batch_size; check CUDA |
| All predictions same value | Constant target | Standardize labels |
| AUC mismatched across folds | Random seed not set | `--seed 42` |
| Test AUC = train AUC | No held-out data | Use scaffold_balanced split |
| Ensemble variance always small | All models identical | Set `--seed` per fold; check randomness |
| SHAP fails on D-MPNN | Custom architecture | Use GradientExplainer not TreeExplainer |
| MolFormer fine-tune slow | All parameters trained | Use LoRA or freeze early layers |
| Calibration broken | Pre-calibrated already | Skip calibration step |

## References

- Yang et al., *J. Chem. Inf. Model.* 59:3370 (2019) -- chemprop D-MPNN.
- Heid et al., *J. Chem. Inf. Model.* 64:9 (2024) -- chemprop 2.0 redesign.
- Wu et al., *Chem. Sci.* 9:513 (2018) -- MoleculeNet benchmark.
- Ross et al. (2022) -- MoLFormer architecture.
- Zhou et al., *Nat. Commun.* 14:3849 (2023) -- Uni-Mol 3D.
- OECD, "OECD Principles for the Validation of QSAR Models" (2007).
- OECD, "(Q)SAR Assessment Framework" (2023).
- Cortés-Ciriano & Bender, *J. Chem. Inf. Model.* 60:1184 (2020) -- conformal prediction for QSAR applicability domain.
- Alperstein et al., *J. Cheminformatics* 15:9 (2023) -- MolBERT chemistry-aware tokenization vs ChemBERTa-2.

## Related Skills

- chemoinformatics/molecular-descriptors - Featurization choices
- chemoinformatics/molecular-standardization - Mandatory upstream
- chemoinformatics/scaffold-analysis - Bemis-Murcko split implementation
- chemoinformatics/admet-prediction - ADMET-specific QSAR
- chemoinformatics/generative-design - QSAR as scoring component
- machine-learning/model-validation - General ML validation principles
- machine-learning/biomarker-discovery - Adjacent ML approaches

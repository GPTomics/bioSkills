---
name: bio-machine-learning-prediction-explanation
description: Explains ML predictions on omics data with SHAP, LIME, and permutation importance, handling the correlated-feature trap, the conditional-vs-interventional Shapley choice, and the attribution-is-not-causation boundary. Use when interpreting an omics classifier, debugging shortcut/batch learning, or deciding whether an attribution ranking can be trusted as biology. For validated feature selection see machine-learning/biomarker-discovery; explanations are not a selection method.
tool_type: python
primary_tool: shap
---

## Version Compatibility

Reference examples tested with: numpy 1.26+, pandas 2.2+, scikit-learn 1.4+, shap 0.44+, lime 0.2+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

Two shap drifts: the Explanation-object plotting API arrived ~0.36 (new-style `shap.plots.*` take an `Explanation`, legacy `shap.summary_plot`/`dependence_plot` take numpy arrays -- mixing them is the most common runtime error); and `TreeExplainer(..., feature_perturbation='auto')` became the default in 0.47 (was `interventional`), so providing or omitting `data=` silently changes the estimand. Always set `feature_perturbation` explicitly. If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt rather than retrying.

# Model Interpretation for Omics Classifiers

**"Which genes drive my classifier?"** -> Compute attributions, but treat them as a description of the model (not biology), choose the Shapley conditioning deliberately, and aggregate over correlated gene modules before ranking.
- Tree models: `shap.TreeExplainer(model, data=background, feature_perturbation='interventional')`
- Model-agnostic local: `lime.lime_tabular.LimeTabularExplainer`
- Model-reliance screen: `sklearn.inspection.permutation_importance`

## The Single Most Important Modern Insight -- Attributions Explain the Model, Not Biology, and Under Correlation the Algorithm Chooses How to Split Credit

A feature attribution describes the function the model learned on this training distribution; it is not a measurement of biology. A high-SHAP gene can be a pure correlate of a batch, scanner, or library-prep signal the model exploited (DeGrave 2021 is the canonical proof). And because genes co-express in tight modules, the attribution algorithm has genuine freedom in how it splits credit among correlated genes -- the choice of conditioning (`tree_path_dependent`/conditional vs `interventional`/marginal) is not a cosmetic knob, it changes *which* genes get credit, and it is a live methodological controversy with no universally correct answer (Janzing 2020). The operational rule: SHAP/LIME rankings are a debugging and hypothesis-generation tool, never a validated biomarker-selection criterion, and within a co-expression module the ordering is not a finding.

## Method Taxonomy

| Method | What it estimates | Correlated-feature behavior | Cost | Best use |
|--------|-------------------|------------------------------|------|----------|
| TreeSHAP `tree_path_dependent` | Conditional Shapley via tree coverage; approximates E[f \| x_S] | Can give nonzero credit to a feature the model never uses (correlation leak); no background needed | Fast, exact for this estimand | Fast cohort summaries when conditional semantics are acceptable |
| TreeSHAP `interventional` | Marginal/do-operator Shapley; features replaced from a background | Zero credit to unused features even if correlated | Scales with background size (~100-1000) | "What the model actually uses"; most defensible default |
| KernelSHAP | Model-agnostic Shapley via masking; assumes independence | Masking lands off-manifold under correlation; corrupted | Expensive | Last resort for non-tree/non-net models |
| DeepSHAP / GradientSHAP | SHAP for nets via backprop relative to a background | Background-dependent under correlation | Moderate | Neural omics models |
| LinearExplainer | Exact Shapley for linear models | `interventional` vs `correlation_dependent` give different values | Cheap | Penalized linear models; choose the mode |
| LIME | Local sparse linear surrogate on perturbed samples | Off-manifold perturbations; unstable across seeds | Moderate | Eyeballing one local prediction, never global ranking |
| Permutation importance | Drop in score when a feature is shuffled | Shuffling A while correlated B intact zeros BOTH; extrapolates | n_repeats x n_features | Global screen on decorrelated features |
| Conditional permutation (Strobl) | Importance permuting within correlated strata | Fairer among correlated predictors | Higher | RF importance under correlation |

## The Correlated-Feature / Conditional-vs-Marginal Problem (the core)

When Shapley values "drop" a feature subset, they replace it by some distribution, and two incompatible choices exist:

- **Conditional / observational** (`tree_path_dependent`): dropped features drawn from `p(x_dropped | x_S)`. A feature the model *never uses* can still receive nonzero attribution purely because it is correlated with a used feature. So **high SHAP does not mean the model relies on this gene.**
- **Marginal / interventional** (`interventional`): dropped features drawn from the marginal `p(x_dropped)`, i.e. `do(x_dropped = background)`. Features the model genuinely ignores get **exactly zero**, even if correlated. Janzing 2020 argues this is the principled "drop" for attribution; it is why modern SHAP added and (pre-0.47) defaulted to it.

There is no free lunch: every method either extrapolates off-manifold (interventional SHAP, unrestricted permutation -- Hooker 2021's "no free variable importance") or leaks credit through correlation (conditional SHAP). Decide which pathology the question can tolerate, and **aggregate attributions over co-expression modules before ranking** -- "gene A ranked above gene B" within a correlated module is governed by off-manifold value-function behavior, not biology.

## Decision Tree by Scenario

| Scenario | Recommended approach | Why |
|----------|---------------------|-----|
| "Which genes is my model actually keying on" (debug shortcuts) | TreeSHAP `interventional` with a representative background | Gives unused genes zero; reveals true reliance |
| "Which genes are informative about the outcome here" (descriptive) | TreeSHAP `tree_path_dependent`, but never call it model reliance | Conditional semantics answer the descriptive question |
| Ranking importance across correlated genes | Aggregate \|SHAP\| within co-expression clusters first | Within-module order is arbitrary |
| Penalized linear model | `LinearExplainer` (choose `interventional` vs `correlation_dependent`) | Or just read the coefficients -- the model is its own explanation |
| One local prediction, communication | LIME or SHAP waterfall, pinned seed, labeled model-internal | Local surrogate; not global, not reproducible across seeds |
| "I want to pick a biomarker panel" | -> machine-learning/biomarker-discovery | SHAP ranking is not validated selection (no FDR, no replication) |
| High-stakes clinical decision | Prefer an inherently interpretable model | Sparse linear/rule list is exact; avoids the conditional-vs-marginal ambiguity (Rudin 2019) |

## SHAP TreeExplainer (set the conditioning explicitly)

**Goal:** Attribute a tree model's predictions to features with a chosen, stated estimand.

**Approach:** Pass a background and `feature_perturbation='interventional'` for "what the model uses," or omit the background and use `tree_path_dependent` for the conditional/descriptive view. The default `'auto'` flips between them based on whether `data=` is given.

```python
import shap
import numpy as np

# Interventional ('what the model uses'): needs a background (~100-1000 rows).
background = shap.utils.sample(X_train, 200)
explainer = shap.TreeExplainer(model, data=background, feature_perturbation='interventional')
sv = explainer(X_test)                                  # modern Explanation object

# Aggregate over correlated modules BEFORE ranking (clusters = a precomputed gene->module map).
mean_abs = np.abs(sv.values).mean(axis=0)
module_importance = {}
for gene, m in zip(X_test.columns, mean_abs):
    module_importance[clusters[gene]] = module_importance.get(clusters[gene], 0) + m
```

## Attribution Is Not Causation, Mechanism, or a Validated Biomarker

Three layers of "not": (1) **not biology** -- the attribution describes the model, which may have exploited a batch/confounder shortcut; (2) **not causation** -- predictive features conflate direct effects, confounders, mediators, and colliders, and turning attribution into a causal claim needs an explicit causal model SHAP does not contain; (3) **not a validated biomarker** -- taking "top-20 SHAP genes" as a panel is same-data feature selection with no FDR and no replication (winner's curse). The strongest *legitimate* use runs the other way: because attribution exposes what the model used, it is one of the best tools to **catch shortcut/batch learning** -- if top features are batch indicators or depth-tracking housekeeping genes, the model is cheating (DeGrave 2021). That is where attribution is most trustworthy.

## LIME and Explanation Instability

LIME fits a sparse linear surrogate to predictions on perturbed samples around one instance. It is **non-reproducible by construction**: different seeds, different `kernel_width`, and `discretize_continuous=True` (the default) each flip the top features, and the per-feature perturbations land off-manifold for correlated genes. Worse, perturbation-based explainers can be deliberately fooled -- a biased model can be wrapped to look innocuous on the out-of-distribution points LIME/KernelSHAP probe (Slack 2020). Use LIME only to eyeball a single prediction's local logic with a pinned seed, never for global ranking, and never as evidence a model is unbiased.

```python
from lime.lime_tabular import LimeTabularExplainer

explainer = LimeTabularExplainer(X_train.values, feature_names=list(X_train.columns),
                                 mode='classification', discretize_continuous=True,
                                 random_state=0)         # pin the seed; still only conditional stability
exp = explainer.explain_instance(X_test.values[0], model.predict_proba, num_features=10, num_samples=5000)
```

## Background / Baseline Choice (the silent attribution-changer)

SHAP explains the deviation from `E[f(X)]` over the background dataset, so the background defines what "absence of a feature" means and changes every attribution. A tumor sample explained against a tumor-heavy vs a healthy-tissue background yields different "important genes" -- only one matches the scientific question. Use a background of real samples representative of the contrast of interest (a single global mean across a heterogeneous cohort is no real sample). `check_additivity=True` (default) raises when the SHAP values plus the base value do not sum to the model output (a local-accuracy violation) -- often a probability-vs-raw-margin or implementation mismatch; investigate rather than disabling it. Attributions in log-odds (`model_output='raw'`) differ from probability space -- report the scale.

## Permutation Importance Also Breaks Under Correlation

A common error is to "fix" SHAP's correlation problem by switching to permutation importance, but it has the same root pathology: shuffling gene A while correlated gene B is intact lets the model recover the signal through B, so both look unimportant (scikit-learn's own multicollinearity example). Unrestricted permutation also forces the model to predict on points that cannot occur (Hooker 2021). For correlated omics predictors, cluster features and keep one per cluster, or use Strobl's conditional permutation (R `party::cforest`, `varimp(conditional=TRUE)`) -- there is no `conditional=` flag in sklearn's `permutation_importance`. Always evaluate permutation importance on held-out data, not training data.

## Per-Method Failure Modes

### Reading high `tree_path_dependent` SHAP as model reliance
- **Trigger:** Using path-dependent SHAP (no background) and concluding the model depends on a top gene.
- **Mechanism:** Conditional Shapley leaks credit to unused-but-correlated features.
- **Symptom:** A gene the model never splits on ranks high.
- **Fix:** Use `interventional` with a background for reliance questions; state the estimand.

### Within-module ranking treated as a finding
- **Trigger:** Reporting "gene A more important than gene B" for co-expressed A, B.
- **Mechanism:** Shapley fairly splits credit, but the split is governed by off-manifold behavior, not biology.
- **Symptom:** The credit split (and sometimes the order) changes with the conditioning mode or when the other module member is dropped.
- **Fix:** Aggregate \|SHAP\| over co-expression modules before ranking.

### SHAP ranking used as feature selection
- **Trigger:** Taking top-k SHAP genes as a biomarker panel.
- **Mechanism:** Same-data selection with no FDR/replication; winner's curse.
- **Symptom:** The panel fails to replicate in an independent cohort.
- **Fix:** Generate hypotheses with SHAP, validate with biomarker-discovery + independent data.

### LIME/KernelSHAP global aggregates trusted
- **Trigger:** Averaging LIME or KernelSHAP across samples for a global ranking.
- **Mechanism:** Seed/kernel instability + off-manifold perturbation; adversarially foolable (Slack 2020).
- **Symptom:** Ranking changes across runs.
- **Fix:** Use TreeSHAP/LinearExplainer for global; keep LIME local.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Background ~100-1000 rows | shap docs | Interventional/Kernel SHAP cost scales with background; shap warns above ~1000 |
| Set `feature_perturbation` explicitly | shap 0.47 changelog | `'auto'` default silently flips estimand on `data=` presence |
| Aggregate over co-expression modules before ranking | Janzing 2020; Aas 2021 | Within-module order is not identifiable from attributions |
| Validate SHAP-derived genes in independent cohorts | Rudin 2019 | Attribution rankings are model-internal, not replicated associations |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `shap.plots.beeswarm` errors on a numpy array | Explanation-object API since ~0.36 | Pass `explainer(X)` (Explanation); use legacy `summary_plot` for arrays |
| Attribution estimand changed silently | `feature_perturbation='auto'` (0.47+) | Set `'interventional'` or `'tree_path_dependent'` explicitly |
| `check_additivity` raises | SHAP values + base do not sum to output (local-accuracy) | Fix the config (raw vs probability); do not just disable |
| `interventional` errors | No `data=` background supplied | Provide a background sample |
| Permutation importance zeros real features | Correlation dilution | Cluster features or use conditional permutation |

## References

- Ribeiro MT, Singh S, Guestrin C. 2016. "Why Should I Trust You?": Explaining the Predictions of Any Classifier. *Proc KDD* 1135-1144.
- Lundberg SM, Lee S-I. 2017. A unified approach to interpreting model predictions. *Adv Neural Inf Process Syst* 30:4765-4774.
- Strobl C, Boulesteix A-L, Kneib T, Augustin T, Zeileis A. 2008. Conditional variable importance for random forests. *BMC Bioinformatics* 9:307.
- Rudin C. 2019. Stop explaining black box machine learning models for high stakes decisions and use interpretable models instead. *Nat Mach Intell* 1:206-215.
- Lundberg SM, Erion G, Chen H, et al. 2020. From local explanations to global understanding with explainable AI for trees. *Nat Mach Intell* 2:56-67.
- Janzing D, Minorics L, Blobaum P. 2020. Feature relevance quantification in explainable AI: a causal problem. *Proc AISTATS* PMLR 108:2907-2916.
- Kumar IE, Venkatasubramanian S, Scheidegger C, Friedler S. 2020. Problems with Shapley-value-based explanations as feature importance measures. *Proc ICML* PMLR 119:5491-5500.
- Slack D, Hilgard S, Jia E, Singh S, Lakkaraju H. 2020. Fooling LIME and SHAP: adversarial attacks on post hoc explanation methods. *Proc AIES*.
- Aas K, Jullum M, Loland A. 2021. Explaining individual predictions when features are dependent: more accurate approximations to Shapley values. *Artif Intell* 298:103502.
- DeGrave AJ, Janizek JD, Lee S-I. 2021. AI for radiographic COVID-19 detection selects shortcuts over signal. *Nat Mach Intell* 3:610-619.
- Hooker G, Mentch L, Zhou S. 2021. Unrestricted permutation forces extrapolation: variable importance requires at least one more model. *Stat Comput* 31:82.

## Related Skills

- machine-learning/omics-classifiers - Train the model being explained; debug batch shortcuts
- machine-learning/biomarker-discovery - Validated feature selection (SHAP ranking is not selection)
- machine-learning/model-validation - Confirm the model generalizes before interpreting it
- data-visualization/heatmaps-clustering - Visualize module-aggregated attributions

---
name: bio-temporal-genomics-temporal-grn
description: Infers directed, time-delayed gene regulatory edges from BULK time-series expression using Granger causality (statsmodels VAR F-test), dynGENIE3 (tree ensembles regressing ODE-derived derivatives; Random Forests by default, Extra-Trees optional), and dynamic Bayesian networks (bnlearn). Use when the output is a RANKED HYPOTHESIS list for perturbation validation, not validated causal edges; deciding Granger vs dynGENIE3 vs DBN by timepoint count and linearity; sizing maxlag against the n>3*maxlag+1 degrees-of-freedom floor; handling stationarity/differencing before Granger; restricting regulators to known TFs; and comparing network rewiring across conditions at matched edge density. Not for single-cell pseudotime GRNs (see gene-regulatory-networks/scenic-regulons) or static co-expression (see gene-regulatory-networks/coexpression-networks).
tool_type: mixed
primary_tool: statsmodels
---

## Version Compatibility

Reference examples tested with: statsmodels 0.14+, numpy 1.26+, pandas 2.2+, dynGENIE3 (GitHub vahuynh/dynGENIE3), bnlearn 4.9+, R 4.x

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Note: bulk time-series GRN inference is low-precision and assumption-heavy. Every edge is a HYPOTHESIS. Results are dominated by the sampling design (interval vs the minutes-scale of transcription, number of timepoints, replicate count), not by the algorithm. A tiny p-value from 6-12 timepoints is not evidence of regulation.

# Temporal Gene Regulatory Network Inference

**"Infer causal regulatory relationships from my time-series expression data"** -> Rank directed, time-delayed TF->target edges from bulk temporal expression, to prioritize perturbation experiments.
- Python: `statsmodels.tsa.stattools.grangercausalitytests()` (VAR F-test on predictive precedence)
- R: `dynGENIE3::dynGENIE3()` (tree ensembles on ODE-derived derivatives); `bnlearn::hc()` + `boot.strength()` (dynamic Bayesian network)

## The governing principle: inference produces ranked HYPOTHESES, not validated causal edges

Bulk temporal GRN inference turns a time course into a ranked list of candidate directed edges whose only honest downstream use is prioritizing perturbation experiments (knockdown / overexpression + re-measure). Two hard facts set the ceiling and must be stated up front, not buried.

1. Granger is PREDICTIVE precedence, not mechanism. It tests whether past X improves prediction of future Y, which is neither necessary nor sufficient for regulation. It collapses in three routine biological situations:
   - Unobserved common driver (confounding). An unmeasured TF, or a shared circadian/cell-cycle oscillation driving hundreds of genes, makes X "Granger-cause" Y with zero direct regulation. Pairwise methods are structurally blind to this; a shared sinusoid manufactures dense, entirely spurious directed structure whose lags are just phase offsets.
   - Sampling coarser than the regulation timescale (aliasing). Transcription acts in minutes; bulk courses are sampled every 1-6 h. When the interval exceeds the regulatory delay, cause and effect land in the same sampled timepoint and directionality becomes unidentifiable. No statistic recovers information the sampling threw away.
   - Non-stationarity. VAR-Granger assumes weak stationarity, but the interesting biology (a stimulus response, a developmental transient, a monotone induction) IS the non-stationary trend, and differencing it away removes the signal (see the differencing dilemma below).
2. Community benchmarks put a LOW ceiling on precision and no single method wins. DREAM5 (Marbach 2012 *Nat Methods* 9:796) evaluated 30+ methods and found time-series network inference is low-precision, no method is best across datasets, and the robust win is the "wisdom of crowds": integrating independent methods beats any one. Prior information (restricting regulators to known TFs) is the other reliable lever.

Operational consequence: restrict regulators to annotated TFs, run more than one method, keep edges recovered by >=2 methods and stable across replicate series, match density before comparing conditions, and hand the top edges to perturbation. This skill is bounded to BULK real-clock-time data; single-cell pseudotime GRN is a different problem (gene-regulatory-networks/scenic-regulons).

## Method selection

| Method | Models | Best when | Fails when |
|--------|--------|-----------|------------|
| Granger (statsmodels) | Bivariate VAR; F-test restricted vs unrestricted | Enough timepoints (n comfortably > 3*maxlag+1); a small a-priori TF->target set; roughly linear, stationary-after-differencing series | 6-12 timepoints (no residual DoF -> no power); genome-wide pairwise (confounding + O(TF*target) tests); saturating/switch-like regulation (linear only) |
| dynGENIE3 (R) | Semi-ODE: trees regress dx/dt on regulator expression | Non-linear / combinatorial regulation; multiple replicates and reasonably dense sampling; a curated regulator list | Sparse or unevenly-spaced timepoints (finite-difference derivative is garbage); calibrated significance is required (it gives a RANKING, no p-values) |
| DBN (bnlearn) | Unrolled first-order Markov Bayesian network across slices | Feedback loops matter (autoregulation, negative feedback); a pre-filtered set of tens-to-low-hundreds of nodes; edge-confidence needed | Genome-wide (super-exponential DAG search); delays longer than one sampling interval (first-order Markov); tiny samples (CI/score tests underpowered) |

Methodology evolves; verify current best practice against each tool's latest documentation before committing to one. The defensible default is to run more than one and intersect.

## Granger causality (Python / statsmodels)

**Goal:** Rank TF->target pairs by whether past TF expression improves prediction of future target expression, with honest multiple-testing control.

**Approach:** Difference all genes uniformly to approach stationarity, select a single lag per pair by BIC (so the reported p-value is not the best-of-several), run ONE F-test at that lag, then BH-correct across pairs. Test only TF->target pairs to shrink the family and encode the TF prior.

The F-test compares an unrestricted VAR (Y on its own lags AND X's lags) to a restricted model (Y on its own lags only); statsmodels reports it as `ssr_ftest`, matching R's `lmtest::grangertest`. Two constraints dominate:

- Degrees-of-freedom floor. After lagging, `n_eff = n - maxlag` rows fit `2*maxlag+1` parameters, so the test is only defined for `n > 3*maxlag + 1`, and barely-defined means no power. With n=8 and maxlag=2 the F-test has ~1 residual DoF: a coin flip. This, not compute, is why genome-wide pairwise Granger fails. Prefer maxlag=1 on short courses.
- Lag selection is itself a multiple test. Taking the minimum p-value over lags 1..maxlag and reporting it as a single test inflates significance. Fix by selecting one lag a priori, or by BIC (below), or by Bonferroni across lags before the across-pairs BH.

```python
import numpy as np
import pandas as pd
from statsmodels.tsa.api import VAR
from statsmodels.tsa.stattools import grangercausalitytests
from statsmodels.stats.multitest import multipletests

# expr_df: genes x timepoints DataFrame; columns MUST be in temporal order.
# Difference uniformly to approach stationarity. Uniform (not per-gene) differencing
# keeps every series on the same footing: mixing I(0) and differenced I(1) series in one
# VAR corrupts the F-test reference distribution. Cost: over-differencing already-stationary
# genes. The deeper tradeoff: differencing removes the trend that CARRIES the regulatory
# signal, so on short courses prefer maxlag=1 over aggressive differencing.
expr_diff = expr_df.diff(axis=1).iloc[:, 1:]

tf_genes = ['TF1', 'TF2', 'TF3']
target_genes = ['geneA', 'geneB', 'geneC']
maxlag = 1  # short courses have ~no DoF beyond lag 1 (need n > 3*maxlag+1)

def granger_pvalue(pair_data, maxlag):
    # column 0 = response Y (target), column 1 = predictor X (TF): tests X -> Y.
    # Select ONE lag by BIC, then run a SINGLE test at it -> avoids the min-p-over-lags
    # multiple test. Guard BIC=0 (no lag structure) up to 1.
    lag = max(1, int(VAR(pair_data).select_order(maxlag).bic))
    res = grangercausalitytests(pair_data, maxlag=[lag])  # list -> tests only this lag
    return res[lag][0]['ssr_ftest'][1], lag  # (p_value, lag)

records = []
for tf in tf_genes:
    for target in target_genes:
        if tf == target:
            continue
        pair = np.column_stack([expr_diff.loc[target].values, expr_diff.loc[tf].values])
        p, lag = granger_pvalue(pair, maxlag)
        records.append({'tf': tf, 'target': target, 'p_value': p, 'lag': lag})

results_df = pd.DataFrame(records)
# multipletests default is Holm-Sidak, NOT BH; force fdr_bh explicitly.
results_df['q_value'] = multipletests(results_df['p_value'], method='fdr_bh')[1]
significant = results_df[results_df['q_value'] < 0.05].sort_values('q_value')
```

Pairwise Granger cannot separate direct regulation from a chain X->Z->Y or a fork Z->{X,Y}. The correct fix is conditional (multivariate) Granger, conditioning on all other regulators' lags, but that explodes the parameter count and is infeasible at transcriptomic sample sizes. Label pairwise output as a CONFOUNDED candidate set, not direct interactions.

## dynGENIE3 (R)

**Goal:** Rank regulator->target edges non-linearly by how much a regulator's current expression predicts a target's temporal derivative.

**Approach:** dynGENIE3 models each gene as `dx_i/dt = f_i(x) - alpha_i * x_i`, estimates `dx_i/dt` by finite differences between consecutive timepoints, and trains a tree ensemble to regress that derivative-plus-decay target on candidate-regulator expression; summed variable importance becomes the edge weight.

```r
library(dynGENIE3)

# TS.data: list of genes x timepoints matrices (one per replicate/series).
# time.points: matching list of time vectors (real deltas -> handles uneven spacing).
expr_list <- list(as.matrix(expr_series1), as.matrix(expr_series2), as.matrix(expr_series3))
time_list <- list(c(0, 4, 8, 12, 24, 48), c(0, 4, 8, 12, 24, 48), c(0, 4, 8, 12, 24, 48))

# Restrict regulators to known TFs (AnimalTFDB / PlantTFDB). This helps TWICE: fewer
# features searched per split (faster) AND a non-TF can never be reported as a regulator
# (higher precision). Single highest-yield precision lever.
tf_indices <- which(rownames(expr_list[[1]]) %in% tf_names)

# tree.method DEFAULTS to 'RF' (Random Forests). Extra-Trees is opt-in: tree.method='ET'
# (the config GENIE3 used to win DREAM4). alpha='from.data' (default) estimates per-gene
# mRNA decay from the data; pass a numeric vector to inject measured half-lives (4sU/BRIC-seq).
res <- dynGENIE3(TS.data = expr_list, time.points = time_list, regulators = tf_indices)

# get.link.list (DOT form) is the dynGENIE3 function. The camelCase getLinkList belongs to
# the separate Bioconductor GENIE3 package -- do not swap them.
link_list <- get.link.list(res$weight.matrix, report.max = 1000)
```

Two properties gate interpretation:
- The weight matrix is a RANKING with no null, no p-value, no calibrated threshold. A "top edge" is top only relative to the others in this run; thresholding by rank (top-K) is unavoidably arbitrary. This is why cross-method agreement and stability matter more here than anywhere.
- Finite-difference derivatives amplify noise. With few, unevenly-spaced timepoints (0,4,8,12,24,48 h is typical) each `dx/dt` rests on one noisy pair and late wide intervals blur short-timescale regulation into a single slope. More REPLICATES (independent derivative samples averaging the noise down) help far more than adding one or two timepoints.

## Dynamic Bayesian networks (R / bnlearn)

**Goal:** Learn a directed network that can represent feedback, with bootstrap edge confidence, over a pre-filtered gene set.

**Approach:** Unroll time into t-1 and t slices and allow edges only from t-1 to t; because A_{t-1}->B_t and B_{t-1}->A_t both point forward, the unrolled graph is acyclic even though the biology has an A<->B feedback loop. So DBNs represent feedback that static Bayesian networks (which must be DAGs) structurally cannot -- the main reason to reach for one. The cost: it is first-order Markov (state at t depends only on t-1; longer delays need t-2/t-3 slices) and the super-exponential DAG search caps realistic inference at tens-to-low-hundreds of nodes, never genome-wide.

```r
library(bnlearn)

# Build the 2-slice frame: columns _t1 (predictors at t-1) and _t (response at t).
n_t <- ncol(expr_mat)
lagged_df <- data.frame(
    t(expr_mat[, 2:n_t]),      # response slice t
    t(expr_mat[, 1:(n_t - 1)]) # predictor slice t-1
)
colnames(lagged_df) <- c(paste0(rownames(expr_mat), '_t'),
                         paste0(rownames(expr_mat), '_t1'))

# Constrain edges to t-1 -> t so the learned graph is a proper DBN transition model.
nodes_t  <- paste0(rownames(expr_mat), '_t')
nodes_t1 <- paste0(rownames(expr_mat), '_t1')
blacklist <- rbind(
    expand.grid(from = nodes_t, to = nodes_t1),   # forbid t -> t-1 (backward in time)
    expand.grid(from = nodes_t1, to = nodes_t1)   # forbid within-slice t-1 edges
)

# score='bic-g': Gaussian BIC; penalizes parameters, guarding the tiny sample against
# overfit. Gaussian assumes linear-Gaussian dependencies (misses threshold logic, like
# Granger); discretizing captures nonlinearity but needs data you do not have on short
# courses. hc is greedy -> trust boot.strength, not one DAG.
boot_res <- boot.strength(lagged_df, R = 200, algorithm = 'hc',
                          algorithm.args = list(score = 'bic-g', blacklist = blacklist))

# strength = fraction of bootstraps containing the arc; direction = fraction of those
# oriented the stated way. direction >= 0.5 is a COIN FLIP -- require >= 0.8 for a
# confidently oriented edge. bnlearn can also compute a data-driven strength threshold:
thr <- attr(boot_res, 'threshold')  # data-driven threshold lives on the bn.strength object, a principled alternative to hand-picked 0.7
confident <- boot_res[boot_res$strength >= max(0.7, thr) & boot_res$direction >= 0.8, ]
```

## Comparing networks across conditions

**Goal:** Identify genuine rewiring between two conditions, not artifacts of threshold choice.

**Approach:** Edge-set differences are dominated by density mismatch and near-threshold flips unless controlled. Compare at MATCHED edge density (top-K from each, same K), and only call an edge gained/lost if it is present-and-bootstrap-stable in one condition and absent-and-stable in the other.

```python
def top_k_edges(edge_df, k):
    return set(map(tuple, edge_df.sort_values('weight', ascending=False)
                   .head(k)[['tf', 'target']].values))

k = min(len(edges_a), len(edges_b))  # density-match BEFORE comparing
set_a, set_b = top_k_edges(edges_a, k), top_k_edges(edges_b, k)
jaccard = len(set_a & set_b) / len(set_a | set_b) if (set_a | set_b) else 0.0
gained, lost = set_b - set_a, set_a - set_b  # keep only bootstrap-stable ones
```

Jaccard heuristics (< 0.3 rewired, > 0.7 conserved) are uncalibrated and, without density-matching, mostly measure the threshold rather than biology -- present them as rough anchors only after matching.

## What experts do instead of trusting one method

- Prior-constrain regulators to annotated TFs (dynGENIE3 `regulators=`; Granger test only TF->target; DBN whitelist/blacklist). Highest-yield, cheapest precision lever.
- Ensemble across methods; edges recovered by >=2 orthogonal methods are the ones worth an experiment (Marbach 2012's wisdom-of-crowds result).
- Require replication across independent time series; bootstrap-subsample and re-rank to separate reproducible edges from artifacts.
- Treat the output as a prioritized hypothesis list for perturbation. Nothing in bulk inference validates an edge; only perturbation does.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| `grangercausalitytests(..., verbose=False)` raises FutureWarning | `verbose` deprecated since statsmodels 0.14, slated for removal | Drop the argument; index the returned dict (`res[lag][0]['ssr_ftest'][1]`) |
| Granger q-values suspiciously optimistic | `min_p` across lags then BH is an uncorrected within-pair multiple test | Fix one lag a priori, or BIC-select one lag then run a single test, or Bonferroni across lags before BH |
| Granger has no power / errors on few timepoints | `n > 3*maxlag+1` barely met -> ~1 residual DoF | Use maxlag=1 on short courses; get more timepoints/replicates before trusting any q-value |
| "dynGENIE3 uses Extra-Trees" | dynGENIE3 defaults to `tree.method='RF'` (Random Forests); ET is opt-in | Pass `tree.method='ET'` if ET is wanted, else describe it as RF |
| dynGENIE3 edges read as calibrated | Importances have no null / no p-value | Threshold by rank explicitly; validate top edges by cross-method agreement + perturbation |
| dynGENIE3 gives garbage on sparse/uneven series | Finite-difference `dx/dt` amplifies noise | Add replicates (independent derivative samples), not just one more timepoint |
| DBN `direction >= 0.5` admits reversed edges | 0.5 = "more often than not" = coin-flip orientation | Require `direction >= 0.8`; consider bnlearn's data-driven strength threshold over a hand-picked 0.7 |
| Pairwise Granger reported as direct regulation | Blind to common drivers / chains; circadian oscillation fabricates dense edges | Label as confounded candidates; restrict to TF->target; intersect methods |
| Jaccard swings wildly between conditions | Density mismatch + near-threshold flips, not biology | Match edge density (top-K each); require bootstrap-stable presence/absence |
| Lag structure vanishes silently | Expression columns not in temporal order | Assert timepoint ordering before any lagging |

## Related Skills

- gene-regulatory-networks/coexpression-networks - Static (non-temporal) co-expression networks
- gene-regulatory-networks/scenic-regulons - Single-cell pseudotime regulon inference (different data and assumptions)
- gene-regulatory-networks/differential-networks - Condition-specific network comparison
- differential-expression/timeseries-de - Filter to temporally-variable genes before edge inference
- data-visualization/network-visualization - Plotting inferred networks

## References

- Granger CWJ. 1969. Investigating causal relations by econometric models and cross-spectral methods. *Econometrica* 37(3):424-438. Predictive-precedence definition of causality.
- Huynh-Thu VA, Geurts P. 2018. dynGENIE3: dynamical GENIE3 for the inference of gene networks from time series expression data. *Sci Rep* 8:3384. Semi-ODE + tree-regression-on-derivative method.
- Huynh-Thu VA, Irrthum A, Wehenkel L, Geurts P. 2010. Inferring regulatory networks from expression data using tree-based methods. *PLoS ONE* 5(9):e12776. GENIE3 tree-based variable selection.
- Marbach D, Costello JC, Kuffner R, et al. 2012. Wisdom of crowds for robust gene network inference. *Nat Methods* 9(8):796-804. Low precision, no single method wins, community-ensemble superiority.
- Scutari M. 2010. Learning Bayesian networks with the bnlearn R package. *J Stat Softw* 35(3):1-22. bnlearn `hc` / `boot.strength` API.
- Friedman N, Murphy K, Russell S. 1998. Learning the structure of dynamic probabilistic networks. *Proc. 14th Conf. on Uncertainty in Artificial Intelligence (UAI)*, pp. 139-147. Score-based DBN structure learning.

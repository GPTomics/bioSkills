# Temporal Gene Regulatory Network Inference - Usage Guide

## Overview

Bulk temporal GRN inference turns a time course into a RANKED LIST OF HYPOTHESES about directed TF->target regulation, not validated causal edges. Granger causality tests predictive precedence (not mechanism) and collapses under unobserved common drivers (an unmeasured TF, or a shared circadian oscillation driving both genes), under sampling coarser than the minutes-scale of transcription, and under the non-stationary transients that are the actual biology of interest. dynGENIE3 gives an importance ranking with no p-values and rests on noise-amplifying finite-difference derivatives. Dynamic Bayesian networks can represent feedback (their real edge over static Bayesian networks) but are first-order Markov and capped at tens of nodes. Community benchmarks show this whole class is low-precision and that no single method wins, so the defensible workflow restricts regulators to known TFs, runs more than one method, keeps edges recovered by multiple methods and stable across replicates, matches density before comparing conditions, and hands the top edges to perturbation experiments.

## Prerequisites

### Python
```bash
pip install statsmodels pandas numpy networkx matplotlib
```

### R
```r
install.packages('bnlearn')
devtools::install_github('vahuynh/dynGENIE3/dynGENIE3R')
```

### Data Requirements
- A time-series expression matrix with columns in temporal order; more timepoints strongly improve lag estimation (Granger power needs n comfortably above 3*maxlag+1).
- A known or predicted transcription factor list (AnimalTFDB, PlantTFDB, or similar) to restrict candidate regulators.
- Multiple biological replicates, which matter more for dynGENIE3 (independent derivative samples) than one extra timepoint.
- Evenly spaced timepoints preferred for Granger; dynGENIE3 handles uneven spacing mechanically but cannot recover dynamics the spacing failed to sample.

## Quick Start

Tell the AI agent what to infer:
- "Infer regulatory relationships between my TFs and targets from time-series expression data"
- "Run Granger causality to find time-delayed gene regulation, with proper multiple-testing"
- "Build a dynamic gene regulatory network from my temporal RNA-seq data"
- "Compare regulatory networks between conditions and tell me which edges are real rewiring"

## Example Prompts

### Granger Causality
> "I have 20 timepoints of RNA-seq data. Test Granger causality between my transcription factors and target genes, selecting the lag by BIC and controlling FDR across pairs."

> "Run pairwise Granger causality on my time-series expression, but only TF->target pairs, and tell me honestly whether my 8 timepoints have enough power."

### dynGENIE3
> "Use dynGENIE3 to infer a gene regulatory network from my time-series expression with a known TF regulator list."

> "I have 3 biological replicates of a developmental time course. Run dynGENIE3 and rank the key regulators."

### Dynamic Bayesian Networks
> "Learn a dynamic Bayesian network from my temporal expression using hill-climbing with BIC, restricted to t-1 to t edges."

> "Bootstrap a dynamic Bayesian network to find confidently oriented regulatory edges among my TFs and targets."

### Network Comparison
> "Compare the regulatory networks between treated and control time courses at matched edge density. Which edges are genuinely gained or lost?"

> "Track network rewiring across my developmental stages and separate real changes from near-threshold flips."

## What the Agent Will Do

1. Load the time-series expression matrix and TF list, and verify the columns are in temporal order.
2. For Granger, difference all genes uniformly to approach stationarity and check the degrees-of-freedom floor.
3. Run the selected inference method: Granger (BIC-selected lag, single test per pair, BH across pairs), dynGENIE3 (RF ensemble on ODE derivatives, regulators restricted to TFs), or a DBN (hill-climbing with a t-1 to t blacklist plus bootstrap edge confidence).
4. Filter edges by q-value, importance rank, or bootstrap strength/direction.
5. Build a directed adjacency matrix and, if requested, compare conditions at matched edge density.
6. Frame the output as a ranked hypothesis list for perturbation, stating the confounding and sampling caveats.

## Tips

- Granger p-values from 6-12 timepoints are essentially untrustworthy: the F-test has almost no residual degrees of freedom. Use maxlag=1 on short courses and treat any hit as a weak prior.
- Do not take the minimum p-value across lags and report it as a single test; that inflates significance. Fix one lag a priori or select it by BIC first.
- Drop the `verbose` argument from `grangercausalitytests`; it is deprecated since statsmodels 0.14.
- Difference uniformly across all genes, never per-gene; mixing differenced and level series corrupts the VAR F-test. Remember differencing removes the trend that carried the regulatory signal.
- dynGENIE3 defaults to Random Forests (`tree.method='RF'`), not Extra-Trees; pass `tree.method='ET'` if Extra-Trees is wanted.
- Restrict regulators to known TFs for both precision and speed; genome-wide pairwise inference is dominated by false positives.
- DBN `direction >= 0.5` is a coin flip; require `direction >= 0.8` for a confidently oriented edge, and consider bnlearn's data-driven strength threshold instead of a hand-picked 0.7.
- Edges recovered by both Granger and dynGENIE3 (or stable across replicate series) are the ones worth an experiment; a single-method top edge is not.
- Compare condition networks at matched edge density (top-K each); raw Jaccard mostly measures the threshold.

## Related Skills

- gene-regulatory-networks/coexpression-networks - Static (non-temporal) co-expression networks
- gene-regulatory-networks/scenic-regulons - Single-cell pseudotime regulon inference (different data and assumptions)
- gene-regulatory-networks/differential-networks - Condition-specific network comparison
- differential-expression/timeseries-de - Filter to temporally-variable genes before edge inference
- data-visualization/network-visualization - Plotting inferred networks

---
name: bio-tcr-bcr-analysis-repertoire-visualization
description: Draws TCR/BCR repertoire figures - V-J chord/circos, CDR3 spectratype, clonal-space stratification, clonal tracking across timepoints, rarefaction/extrapolation curves, overlap heatmaps, and clonotype-similarity networks - and encodes how to read them. Use when choosing between a raw Shannon bar and a rarefaction curve for a diversity comparison; deciding a depth-robust overlap metric (Morisita-Horn) vs a set metric (Jaccard) for a heatmap; setting the distance threshold that defines a clonotype-similarity network; interpreting a Gaussian vs skewed spectratype as polyclonal vs clonally expanded; or laying out clonal-space and clone-tracking plots. Covers VDJtools PlotFancyVJUsage/RarefactionPlot, R circlize and iNEXT, and matplotlib/seaborn recipes.
tool_type: mixed
primary_tool: VDJtools
goal_approach_exempt: true
---

## Version Compatibility

Reference examples tested with: matplotlib 3.8+, seaborn 0.13+, pandas 2.2+, numpy 1.26+; R circlize 0.4+, iNEXT 3.0+; VDJtools 1.2.1+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- R: `packageVersion('<pkg>')` then `?function_name` to verify parameters
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Note: matplotlib 3.9 removed `plt.cm.get_cmap(name, N)`; use `plt.get_cmap(name)` and sample it, or `matplotlib.colormaps[name].resampled(N)`. seaborn 0.13 deprecated bare `palette=` without `hue=`; assign `hue=` and `legend=False`.

# Repertoire Visualization

**"Plot my TCR/BCR repertoire"** -> Render V-J usage, CDR3-length spectratype, clonal-space, clonal tracking, diversity/rarefaction, overlap, and similarity-network figures, and read each one correctly.
- CLI: `vdjtools PlotFancyVJUsage`, `vdjtools RarefactionPlot` (delegates plotting to R)
- R: `circlize::chordDiagram` (V-J), `iNEXT::iNEXT` + `ggiNEXT` (rarefaction/extrapolation)
- Python: `matplotlib`/`seaborn` for bespoke figures; `numpy` multinomial resampling for rarefaction; `networkx` + `rapidfuzz` for similarity networks

## The governing principle: every figure inherits two upstream choices

A repertoire figure is only as valid as the numbers behind it, and those numbers depend on two decisions made before any plot is drawn:

1. The clonotype definition. A clonotype = CDR3 (nucleotide OR amino acid) + optionally V and J. Amino-acid CDR3 collapses convergent recombination (many nt rearrangements -> one aa clonotype), inflating apparent sharing and deflating richness; nt CDR3 is more conservative (Venturi et al. 2006 *PNAS* 103:18691). Every count on every axis shifts with this choice, so it must be stated in the figure and held constant across all samples in a comparison.
2. The sequencing depth. Richness, Shannon, clonality, and set-overlap (Jaccard, shared-clonotype counts) are strictly increasing functions of reads sequenced - the rare-clone tail never saturates. Diversity, rarefaction, and overlap figures are comparable across samples ONLY after depth normalization: downsample all samples to a common depth, or read rarefaction curves at a shared x. Comparing raw values across libraries of unequal depth measures sequencing effort, not biology (Chao et al. 2014 *Ecol Monogr* 84:45; Greiff et al. 2015 *Genome Med* 7:49). This caveat governs the diversity, rarefaction, and overlap recipes below.

## Choosing the figure and reading it

| Figure | What it reveals | How to compare correctly |
|--------|-----------------|--------------------------|
| V-J chord/circos | Combinatorial V-J pairing bias within one sample | Descriptive per sample; never compare raw usage across primer sets/platforms (primer bias masquerades as biology) |
| CDR3 spectratype | Clonal structure: Gaussian length distribution = polyclonal/naive, spikes = expansion | Weight by frequency to see expansions; by unique clonotypes to see underlying diversity |
| Clonal-space / proportion | Fraction of repertoire held by rare vs expanded clones | Bin by clone frequency (strata), not raw count; robust to depth if frequency-based |
| Clonal tracking | Expansion/contraction/persistence of clones over time | Downsample timepoints to common depth first; an "absent" clone is often a sampling zero |
| Rarefaction/extrapolation | Diversity comparison done right (interpolate to common depth) | Read all curves at a shared x; extrapolate at most 2-3x observed depth |
| Overlap heatmap | Pairwise repertoire similarity | Use Morisita-Horn (depth-robust) on depth-normalized samples; Jaccard is depth-biased |
| Similarity network | Clusters of related CDR3s (candidate specificity groups) | Structure depends entirely on the distance threshold; state it and test sensitivity |

## V-J usage: chord/circos and heatmap

The chord/circos shows which V and J segments pair, weighted by clonotype count, within one sample.

VDJtools route (delegates plotting to R; requires `RInstall` once):

```bash
vdjtools PlotFancyVJUsage -m metadata.txt output_dir/
```

R with circlize - a V-by-J count matrix becomes a chord diagram:

```r
library(circlize)

plot_vj_chord <- function(clone_df) {
    vj_matrix <- table(clone_df$v_gene, clone_df$j_gene)
    chordDiagram(vj_matrix, transparency = 0.5, annotationTrack = c('grid', 'name'))
}
```

Python heatmap alternative (samples x V gene, or V x J for one sample) avoids a chord dependency and reads more quantitatively for many segments:

```python
import seaborn as sns
import matplotlib.pyplot as plt

def plot_vj_heatmap(clone_df):
    vj = clone_df.pivot_table(index='v_gene', columns='j_gene', values='frequency', aggfunc='sum', fill_value=0)
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.heatmap(vj, cmap='viridis', ax=ax)
    ax.set_title('V-J pairing frequency')
    return fig
```

## CDR3 spectratype: length distribution reveals clonal structure

The spectratype is a histogram of CDR3 length. A roughly Gaussian, bell-shaped distribution indicates a diverse, polyclonal (naive-like) repertoire; skew or discrete spikes at particular lengths indicate clonal expansion(s). The classic immunoscope spectratype (and VDJtools `CalcSpectratype`) bins nucleotide length, where in-frame clones sit 3 nt apart and the periodicity is part of the readout; amino-acid length is a common simplification - state which is plotted. Weighting by read/UMI frequency shows expansions; weighting by unique clonotypes shows the underlying diversity - the two views can look opposite, so label which is plotted.

```python
def plot_spectratype(clone_df, length_col='cdr3_length'):
    fig, ax = plt.subplots(figsize=(9, 5))
    bins = range(clone_df[length_col].min(), clone_df[length_col].max() + 2)
    ax.hist(clone_df[length_col], bins=bins, weights=clone_df['frequency'], color='steelblue')
    ax.set_xlabel('CDR3 length (aa)')
    ax.set_ylabel('Frequency')  # Gaussian = polyclonal; spikes = clonal expansion
    return fig
```

## Clonal-space / proportion: how much repertoire the big clones hold

Bin clones into frequency strata (Rare / Small / Medium / Large / Hyperexpanded) and stack their summed frequency. Because strata are defined on frequency, this view is comparatively depth-robust and shows clonal-space homeostasis at a glance. A treemap of top clones is an alternative when individual dominant clones matter.

```python
import numpy as np

def plot_clonal_space(clone_df, sample_col='sample'):
    # Frequency-based strata (immunarch homeo convention); frequency makes this depth-robust
    edges = [0, 1e-4, 1e-3, 1e-2, 1e-1, 1.0]
    labels = ['Rare', 'Small', 'Medium', 'Large', 'Hyperexpanded']
    clone_df = clone_df.copy()
    clone_df['stratum'] = pd.cut(clone_df['frequency'], bins=edges, labels=labels)
    space = clone_df.groupby([sample_col, 'stratum'], observed=True)['frequency'].sum().unstack(fill_value=0)
    fig, ax = plt.subplots(figsize=(8, 5))
    space[labels].plot(kind='bar', stacked=True, ax=ax, colormap='viridis')
    ax.set_ylabel('Fraction of repertoire')
    return fig
```

## Clonal tracking across timepoints

Track individual clone frequencies over an ordered sample set (vaccination, infection, therapy). A line plot of the top clones, or an alluvial/stream for the same data, shows expansion, contraction, and persistence. Downsample timepoints to a common depth before declaring contraction - a clone scored absent at one timepoint is frequently below the sampling floor, not truly gone.

```python
def plot_clone_tracking(clone_df, top_n=10, clone_col='cdr3_aa', time_col='timepoint'):
    top = clone_df.groupby(clone_col)['frequency'].sum().nlargest(top_n).index
    fig, ax = plt.subplots(figsize=(9, 5))
    for clone in top:
        d = clone_df[clone_df[clone_col] == clone].sort_values(time_col)
        ax.plot(d[time_col], d['frequency'], marker='o', label=clone[:12])
    ax.set_xlabel('Timepoint'); ax.set_ylabel('Clone frequency')
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=7)
    return fig
```

## Rarefaction/extrapolation: the correct way to compare diversity

A bare bar of Shannon (or observed richness) across samples of unequal depth is misleading - it plots depth as much as biology. The defensible comparison is a rarefaction/extrapolation curve: interpolate each sample down to (and extrapolate modestly above) a shared depth, then read diversity at a common x. iNEXT computes this for Hill numbers q=0/1/2 with confidence intervals (Hsieh et al. 2016 *Methods Ecol Evol* 7:1451; Chao et al. 2014 *Ecol Monogr* 84:45).

R route (gold standard, gives CIs):

```r
library(iNEXT)

# count_list: named list of per-sample integer clonotype-count vectors
out <- iNEXT(count_list, q = c(0, 1, 2), datatype = 'abundance')
ggiNEXT(out, type = 1)  # diversity vs sample size, read at a common x
```

VDJtools route: `vdjtools RarefactionPlot -m metadata.txt output_dir/`.

Python resampling route (interpolation by multinomial subsampling) when iNEXT is unavailable:

```python
def rarefaction_curve(counts, depths, reps=20, rng=None):
    # Interpolate observed richness by drawing 'm' reads without-replacement-like via multinomial
    rng = rng or np.random.default_rng(0)
    counts = np.asarray(counts, dtype=float)
    p = counts / counts.sum()
    total = int(counts.sum())
    richness = []
    for m in depths:
        if m > total:  # interpolation only; do not extrapolate past observed depth here
            richness.append(np.nan); continue
        obs = [np.count_nonzero(rng.multinomial(m, p)) for _ in range(reps)]
        richness.append(np.mean(obs))
    return richness
```

## Overlap heatmap: pick a depth-robust metric

An N-by-N similarity heatmap summarizes pairwise repertoire overlap. The metric choice is the decision: Morisita-Horn is abundance-weighted and near-invariant to sample size, so it is the default across unequal depths; Jaccard (presence/absence) is dominated by the shallower sample's depth and should not be compared across unequal-depth pairs. Compute the matrix in vdjtools-analysis on depth-normalized samples, then render it here. Label the metric in the title.

```python
def plot_overlap_heatmap(overlap_matrix, metric='Morisita-Horn'):
    fig, ax = plt.subplots(figsize=(7, 6))
    sns.heatmap(overlap_matrix, annot=True, fmt='.2f', cmap='YlOrRd', vmin=0, vmax=1, square=True, ax=ax)
    ax.set_title(f'Repertoire overlap ({metric})')  # state the metric; Jaccard is depth-biased
    return fig
```

## Clonotype-similarity network

Nodes are CDR3s, edges connect sequences within a chosen distance, and clusters are candidate specificity groups. The network structure depends entirely on the threshold: too loose chains distinct clones together, too tight fragments real groups. State the threshold and metric, and test sensitivity before interpreting clusters. Same-length CDR3s with Hamming distance is the conservative default; Levenshtein allows indels.

```python
import networkx as nx
from rapidfuzz.distance import Levenshtein

def build_similarity_network(clone_df, max_norm_dist=0.15, clone_col='cdr3_aa'):
    # normalized_similarity in [0,1]; edge when 1 - similarity <= threshold. Structure is threshold-dependent.
    clones = clone_df[clone_col].unique()
    g = nx.Graph()
    g.add_nodes_from(clones)
    for i, a in enumerate(clones):
        for b in clones[i + 1:]:
            if 1 - Levenshtein.normalized_similarity(a, b) <= max_norm_dist:
                g.add_edge(a, b)
    return g
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Diversity "differs" between groups but tracks library size | Bare Shannon/richness bar across unequal depth | Downsample to common depth or plot rarefaction curves read at a shared x (iNEXT/RarefactionPlot) |
| Overlap heatmap shows huge differences driven by one shallow sample | Jaccard/shared-count on raw, unequal-depth counts | Use Morisita-Horn on depth-normalized samples; label the metric |
| Network clusters change completely on re-run with a new cutoff | Arbitrary distance threshold; single-linkage chaining | Fix and report the metric + threshold; run a sensitivity sweep; prefer Hamming on same-length CDR3s |
| Two figures disagree on sharing/diversity | Built on different clonotype definitions (nt vs aa, +/- V/J) | Hold one clonotype definition constant across all figures in a comparison and state it |
| Spectratype "expansion" vanishes when re-plotted | Switched between frequency-weighted and clonotype-weighted histogram | Choose one weighting per figure and label it (frequency shows expansions) |
| V-usage differences between batches look biological | Multiplex-PCR primer bias | Compare usage only within one protocol, or use UMI/5'RACE data |
| `plt.cm.get_cmap(name, N)` raises AttributeError | Removed in matplotlib 3.9+ | Use `plt.get_cmap(name)` or `matplotlib.colormaps[name].resampled(N)` |

## Related Skills

- vdjtools-analysis - Compute the diversity/overlap inputs (depth-normalized)
- mixcr-analysis - Produce clonotype tables
- scirpy-analysis - Single-cell clonal-expansion overlays
- data-visualization/ggplot2-fundamentals - General ggplot2 grammar
- data-visualization/heatmaps-clustering - Overlap/usage heatmap techniques

## References

- Shugay M, et al. VDJtools: unifying post-analysis of T cell receptor repertoires. *PLoS Comput Biol* 2015; 11(11):e1004503. (PlotFancyVJUsage, RarefactionPlot, overlap metrics.)
- Chao A, Gotelli NJ, Hsieh TC, Sander EL, Ma KH, Colwell RK, Ellison AM. Rarefaction and extrapolation with Hill numbers. *Ecol Monogr* 2014; 84(1):45-67. (Interpolation/extrapolation framework.)
- Hsieh TC, Ma KH, Chao A. iNEXT: an R package for rarefaction and extrapolation of species diversity (Hill numbers). *Methods Ecol Evol* 2016; 7(12):1451-1456. (Rarefaction/extrapolation curves with CIs.)
- Greiff V, et al. A bioinformatic framework for immune repertoire diversity profiling. *Genome Med* 2015; 7:49. (Hill-number diversity profiling of repertoires.)
- Venturi V, et al. Sharing of T cell receptors in antigen-specific responses is driven by convergent recombination. *PNAS* 2006; 103(49):18691-18696. (aa-clonotype sharing inflated by convergence.)
- Chao A. Nonparametric estimation of the number of classes in a population. *Scand J Stat* 1984; 11:265-270. (Chao1 richness estimator.)

---
name: bio-tcr-bcr-analysis-vdjtools-analysis
description: Computes immune-repertoire diversity, clonal structure, overlap, and segment usage from TCR/BCR clonotype tables with VDJtools (immunarch as the modern R alternative). Use when deciding which diversity estimator answers a question (q=0 observed richness/chao1/chaoE, q=1 shannonWienerIndex, q=2 inverseSimpson as a Hill profile); normalizing sequencing depth before any cross-sample claim (DownSample or the resampled CalcDiversityStats table); choosing an overlap metric (depth-robust MorisitaHorn/F2 vs depth-biased Jaccard/public counts) and a clonotype match key (-i nt/aa, +/-V/J); summarizing clonality as 1 - normalizedShannonWienerIndex; reading spectratype and V-J usage under primer bias; interpreting public clonotypes; and choosing VDJtools (stable Java CLI) vs immunarch (active tidy R).
tool_type: cli
primary_tool: VDJtools
---

## Version Compatibility

Reference examples tested with: VDJtools 1.2.1+, Java (JRE 8+), R 4.x with ggplot2/reshape2/gridExtra (for Plot* modules), immunarch 0.9+/1.0+

Before using code patterns, verify installed versions match. If versions differ:
- CLI: `java -jar vdjtools.jar` prints the current routine list; `<routine>` with no args prints its flags
- R: `packageVersion('immunarch')` then `?repDiversity` / `?repOverlap` to confirm `.method` strings

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Note: routine names are CamelCase and case-sensitive, and the depth-resampling routine is `DownSample` (capital S). Run `RInstall` once so the `Plot*` modules can call R. VDJtools is post-analysis only: it consumes clonotype tables (from MiXCR etc.), not FASTQ.

# VDJtools Analysis

**"Compute diversity and compare my TCR/BCR repertoires"** -> summarize each repertoire's clonal structure, compare samples at equal depth, and quantify overlap.
- CLI: `java -jar vdjtools.jar CalcDiversityStats | CalcPairwiseDistances | TrackClonotypes`
- R alternative: immunarch `repDiversity()`, `repOverlap()`, `repClonality()`

## The governing principle: diversity is sampling-depth-dependent

A repertoire is a sample of an enormous, unevenly expanded clonal population with a long tail of rare clonotypes, so observed richness never saturates: deeper sequencing keeps discovering new clonotypes. Observed richness, Shannon entropy, clonality, Jaccard, and shared-clonotype counts are all functions of read depth. Comparing raw values across libraries of unequal depth measures depth, not biology -- this is the field's single most common and most invalidating error.

The fix is mandatory before any cross-sample claim: bring all samples to a common depth. Two routes:
- `DownSample -x <reads>` every sample to a shared depth, then analyze; or
- read `CalcDiversityStats` at a common depth from its resampled table. `CalcDiversityStats` emits two tables, `diversity.<i>.txt` (original) and `diversity.<i>.resampled.txt` (downsampled to the smallest sample or `-x`). Use the resampled/normalized values for between-sample comparison; the original table is for within-sample description only.

Choosing the normalization depth is itself a decision: downsampling every sample to the cohort minimum discards data and can leave everyone underpowered if one library is tiny. Set the common depth near the cohort's lower quartile, and EXCLUDE (do not drag everyone down to) any sample far below it -- a sample whose rarefaction curve is still steeply climbing well below the chosen depth is under-sampled and cannot support a diversity claim at all. Report the chosen depth and any excluded samples. `PlotQuantileStats` and the rarefaction curves show which samples are safe to include.

Rarefaction makes the problem visible: `RarefactionPlot` draws interpolated + extrapolated diversity-vs-depth curves. Compare samples at a common x, never at curve endpoints of different depth. Extrapolation is reliable only to ~2-3x observed depth and degrades for q=0 (Chao 2014).

## Report a Hill profile, not one number

A single index misleads because indices weight the abundance distribution differently. Report the Hill profile -- effective number of clonotypes at orders q=0, 1, 2 -- whose shape (steep drop from q=0 to q=2 = a few dominant clones over a large rare tail) is the informative object (Greiff 2015 *Genome Med* 7:49; Chao 2014 *Ecol Monogr* 84:45-67). Two repertoires can share richness yet have opposite clonality.

`CalcDiversityStats` emits these columns (each with `_mean`/`_std`); Gini is NOT among them (it is an immunarch option, not a VDJtools output):

| Column | Hill order | Question it answers | Depth-robustness |
|--------|-----------|---------------------|------------------|
| observedDiversity | q=0 | How many distinct clonotypes were seen | Poor -- must downsample |
| chao1 | q=0 | Nonparametric richness lower bound (uses singletons f1, doubletons f2) | Poor; breaks without count data (f2=0) or if rare clones were pre-filtered; PCR error inflates it |
| chaoE | q=0 | Chao richness extrapolated, normalized for cross-sample use | Moderate (VDJtools' preferred richness proxy) |
| efronThisted | q=0 | Efron-Thisted lower-bound total diversity | Poor; a lower bound, not the truth |
| shannonWienerIndex | q=1 | exp(Shannon), effective number weighting by frequency | Moderate |
| normalizedShannonWienerIndex | -- | Pielou evenness H'/ln(S), range 0-1 | Depth-dependent through ln(S) |
| inverseSimpson | q=2 | 1/sum(p^2), dominated by abundant clones | Best -- most depth-robust |
| d50 | -- | Fewest top clones covering 50% of reads | Poor; coarse descriptor |

Default: report q=0 (chaoE or downsampled observedDiversity), q=1 (shannonWienerIndex), and q=2 (inverseSimpson) together. Feed chao1/efronThisted only genuine count data with singletons and doubletons; on non-UMI, non-error-corrected data, PCR/sequencing errors manufacture singletons and inflate them arbitrarily.

## Clonality: the field default and its three flaws

Clonality = 1 - normalizedShannonWienerIndex = 1 - H'/ln(S). It runs 0 (even/polyclonal) to 1 (one clone dominates) and is the near-universal one-number summary because it is bounded and intuitive. State its flaws in any report:
1. Depth-dependent through the ln(S) denominator -- compare clonality only on depth-normalized samples.
2. It discards richness (it is a rescaled evenness): two repertoires with identical clonality can differ 100x in richness.
3. It is dominated by the middle of the abundance distribution, not the top clones a clinician cares about.

Fix: report clonality alongside a q=2 Hill number (inverseSimpson) and a rarefaction curve, never alone.

## Overlap: pick a depth-robust metric and hold the match key fixed

`CalcPairwiseDistances` builds an N x N matrix; `OverlapPair` compares two samples. Any count-of-shared-clonotypes or set index is dominated by the shallower sample's depth: a clone can only be shared if sampled in both, so the shallow sample caps the intersection. Downsample both samples to a common depth first, and prefer abundance-weighted metrics.

| Metric | Basis | Best when | Fails when |
|--------|-------|-----------|------------|
| MorisitaHorn | Abundance, size-normalized | Unequal depth; the default choice | -- (near-invariant to depth; dominated by abundant shared clones) |
| F2 | Sum of per-clonotype geometric-mean frequencies | Frequency-weighted overlap robust to a single dominant shared clone | -- (preferred VDJtools frequency metric) |
| F | Geometric mean of summed shared frequencies | Quick frequency overlap | One large shared clone dominates it |
| R | Pearson of log-frequencies over shared clones only | Concordance of abundances among shared clones | Ignores private clones entirely |
| D | Shared count / geometric-mean diversities | Descriptive | Numerator (shared count) still depth-biased |
| Jaccard | Presence/absence | Equal-depth, denoised samples only | Dominated by the shallower sample's depth |

Overlap magnitude also swings by orders of magnitude with the clonotype match key (`-i`): `nt` (strict, few coincidental shares) vs `aa` (convergent recombination inflates sharing), and whether V/J must match (`ntV`, `ntVJ`, `aaVJ`, ...). Fix one key and hold it constant across every comparison in a study; state it in every figure. `TrackClonotypes` does ordered all-vs-all intersection for time courses -- a clone scoring "absent" at a timepoint is often a sampling zero, so downsample timepoints to common depth before declaring contraction.

An overlap number is only interpretable against a null: some sharing is expected by chance from convergent recombination of high-Pgen clonotypes. To claim overlap EXCEEDS chance, compare the observed statistic to a background of unrelated-donor pairs, or to shuffled/label-permuted repertoires at the same depth, and for public-clonotype claims condition on generation probability (specificity-annotation). Two related individuals or two timepoints from one host will always overlap more than two random donors regardless of biology.

## Public is not antigen-driven

"Public" (a clonotype shared across individuals) is largely an artifact of generation probability, not shared antigen selection. High-Pgen CDR3s -- short, few insertions, near-germline (and fetal-generated) -- are independently produced by many donors, and convergent recombination compounds this at the aa level (Venturi 2006 *PNAS* 103:18691). So a public/shared count is enriched for stochastic high-Pgen sequences, not evidence of a shared response. To argue antigen association, condition on Pgen (OLGA/IGoR) or intersect with an antigen database (`ScanDatabase` against VDJdb; immunarch `dbAnnotate` against VDJdb/McPAS-TCR) -- and even a database hit is a sequence match, not proof of binding. Hand Pgen-aware interpretation off to specificity-annotation.

## Segment usage and spectratype

`CalcSegmentUsage` yields per-sample V/J frequency vectors; `CalcSpectratype` yields the CDR3-length histogram (`PlotFancySpectratype` overlays the top-N clones; `PlotSpectratypeV` stacks by V family; `PlotFancyVJUsage` is the V-J chord plot).
- Spectratype shape: a Gaussian/bell length distribution indicates a diverse polyclonal (naive-like) repertoire; skew or spikes at particular lengths indicate clonal expansion(s). Weighting by reads shows expansions; weighting by unique clonotypes shows underlying diversity.
- Confound: multiplex-PCR primer sets have V-gene-specific amplification bias, so apparent V/J usage differences between platforms or batches are frequently primer artifacts, not biology (Barennes 2021 *Nat Biotechnol* 39:236). Compare usage only within one protocol, or use 5'-RACE/UMI data. Usage vectors are compositional (sum to 1): CLR-transform before PCA and check that PC1 is not just depth/batch.

## VDJtools vs immunarch

Both consume the same clonotype tables; pick by pipeline, not by metric.

| | VDJtools | immunarch |
|--|----------|-----------|
| Language | Java CLI (calls R for plots) | R / tidyverse |
| Maintenance | Stable, low activity (~1.2.1) | Actively maintained (v1.0 adds `airr_*`) |
| Ingestion | `Convert -S <fmt>` | `repLoad()` auto-detects MiXCR/Adaptive/10x/AIRR/VDJtools |
| Plotting | Fixed `Plot*` PDFs | `vis()` returns editable ggplot objects |
| Strengths | Reference F/F2/chaoE + resampled tables; legacy reproducibility; pairs with MiXCR/VDJdb | 10x single-cell, k-mer/motif, publication plots, ML feature matrices |

Prefer VDJtools for CLI/legacy MiXCR pipelines and its exact resampled diversity tables; prefer immunarch for R, single-cell, k-mer/motif, or editable figures. Many groups convert with VDJtools and analyze/plot with immunarch.

Core operations in immunarch (verify `.method` strings on the installed version):

```r
library(immunarch)
data <- repLoad('samples_dir/')                       # metadata + tidy clonotype tables

repDiversity(data$data, .method = 'raref')            # rarefaction/extrapolation curves (the depth control)
repDiversity(data$data, .method = 'hill')             # Hill profile across q
repDiversity(data$data, .method = 'inv.simp')         # q=2, depth-robust
repClonality(data$data, .method = 'homeo')            # clonal-space homeostasis (Rare..Hyperexpanded bins)
repOverlap(data$data, .method = 'morisita')           # depth-robust overlap; 'jaccard'/'public' are depth-biased
geneUsage(data$data[[1]])                             # V/J usage vector
trackClonotypes(data$data, list('Sample1', 1:10))     # longitudinal tracking
dbAnnotate(data$data, vdjdb, 'CDR3.aa', 'cdr3')       # antigen-database annotation
```

## Prepare and normalize (CLI)

Convert upstream output, drop nonfunctional clones, and downsample to a shared depth before any comparison.

```bash
# Import MiXCR clonotypes to VDJtools format (also: -S migec/immunoseq/imgt/vidjil ...)
java -jar vdjtools.jar Convert -S mixcr mixcr_clones.txt converted/

# Keep only functional (in-frame, no-stop) clonotypes for functional-repertoire analysis
java -jar vdjtools.jar FilterNonFunctional -m metadata.txt filtered/

# Optional: remove cross-sample contamination (barcode switching / chimeras) before cross-sample work
java -jar vdjtools.jar Decontaminate -m metadata.txt decontaminated/

# Downsample every sample to a common read depth (capital S; -x/--size = target reads)
java -jar vdjtools.jar DownSample -x 100000 -m metadata.txt downsampled/
```

The metadata file is tab-delimited: a `#file.name  sample.id  <covariate...>` header row, then one row per sample. Most multi-sample routines consume it via `-m`.

## Diversity and overlap (CLI)

```bash
# Diversity: emits diversity.<i>.txt (original) AND diversity.<i>.resampled.txt (depth-normalized)
# Compare across samples using the RESAMPLED table only.
java -jar vdjtools.jar CalcDiversityStats -m metadata.txt diversity/

# Rarefaction curves -- the correct visual for comparing diversity across depths (read at common x)
java -jar vdjtools.jar RarefactionPlot -m metadata.txt rarefaction/

# Pairwise overlap; -i sets the clonotype match key (hold constant study-wide).
# Report MorisitaHorn / F2 columns; treat Jaccard as depth-biased.
java -jar vdjtools.jar CalcPairwiseDistances -i aa -m metadata.txt overlap/
java -jar vdjtools.jar ClusterSamples -e MorisitaHorn overlap/ clustered/

# Longitudinal tracking across an ordered sample set (downsample timepoints first)
java -jar vdjtools.jar TrackClonotypes -m metadata_timecourse.txt tracking/
```

## Parse VDJtools output in Python

```python
import pandas as pd

def load_resampled_diversity(prefix):
    '''Load the depth-normalized diversity table for cross-sample comparison.'''
    return pd.read_csv(f'{prefix}.strict.resampled.txt', sep='\t')

def load_overlap_matrix(path, metric='MorisitaHorn'):
    '''Load one depth-robust overlap metric from the pairwise-distance output.'''
    df = pd.read_csv(path, sep='\t')
    return df.pivot(index='1_sample_id', columns='2_sample_id', values=metric)
```

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Deeper libraries look 'more diverse' every time | Comparing raw richness/Shannon/clonality across unequal depth -- measuring depth, not biology | `DownSample` to a common depth or read the `.resampled.txt` table; compare rarefaction curves at a common x |
| Overlap flips when samples are swapped or re-sequenced | Jaccard / public counts dominated by the shallower sample's depth | Downsample both, and report MorisitaHorn or F2 |
| Overlap magnitude differs wildly between studies | Different clonotype match key (`-i` aa vs nt, +/-V/J) | Fix one `-i` value and state it in every figure |
| chao1/efronThisted are NaN or absurdly large | Fed frequency-only or rare-clone-filtered data, or PCR errors created singletons | Provide genuine count data (singletons/doubletons); UMI/error-correct first; or use inverseSimpson |
| One clonality number reported as 'the diversity' | Clonality is a rescaled evenness that discards richness and is depth-dependent | Report Hill q=0/1/2 (add inverseSimpson) plus a rarefaction curve |
| 'Public' clones claimed as antigen-specific | Publicity is mostly high-Pgen convergent recombination | Condition on Pgen (OLGA) or intersect VDJdb via `ScanDatabase`; hand off to specificity-annotation |
| V/J usage differs between cohorts by platform | Multiplex-PCR primer bias, not biology | Compare usage only within a protocol; CLR-transform before PCA and check PC1 is not depth/batch |
| `Plot*` routine errors on start | R plotting dependencies missing | Run `java -jar vdjtools.jar RInstall` once |

## Related Skills

- mixcr-analysis - Generate input clonotype tables
- repertoire-visualization - Rarefaction, spectratype and overlap figures
- immcantation-analysis - BCR-aware diversity and clonal analysis
- specificity-annotation - Pgen-aware interpretation of public clonotypes
- experimental-design/sample-size - Sequencing depth and power planning
- workflows/tcr-pipeline - End-to-end orchestration

## References

- Shugay M, et al. VDJtools: unifying post-analysis of T cell receptor repertoires. *PLoS Comput Biol* 2015; 11(11):e1004503.
- Chao A, et al. Rarefaction and extrapolation with Hill numbers: a framework for sampling and estimation in species diversity studies. *Ecol Monogr* 2014; 84(1):45-67.
- Greiff V, et al. A bioinformatic framework for immune repertoire diversity profiling enables detection of immunological status. *Genome Med* 2015; 7:49.
- Chao A. Nonparametric estimation of the number of classes in a population. *Scand J Stat* 1984; 11:265-270.
- Venturi V, et al. Sharing of T cell receptors in antigen-specific responses is driven by convergent recombination. *PNAS* 2006; 103(49):18691-18696.
- Barennes P, et al. Benchmarking of T cell receptor repertoire profiling methods reveals large systematic biases. *Nat Biotechnol* 2021; 39:236-245.
- ImmunoMind Team. immunarch: an R package for painless analysis of T-cell and B-cell immune repertoires. CRAN / immunarch.com (v1.x).

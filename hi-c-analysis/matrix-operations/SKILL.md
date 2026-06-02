---
name: bio-hi-c-analysis-matrix-operations
description: Balances Hi-C contact matrices (ICE via cooler.balance_cooler, KR/SCALE/VC context), computes distance-decay expected with cooltools (expected_cis per-diagonal P(s), expected_trans scalar), builds observed/expected (O/E) matrices, and diagnoses polymer state from the P(s) log-derivative. Covers the within-matrix-vs-cross-sample distinction (balancing is NOT a normalizer), the equal-visibility assumption that CNV/aneuploidy violates (use raw counts for copy-number), cis-only balancing, mad_max/blacklist masking before balancing, multiplicative cooler weights vs divisive juicer weights, and the resolution-vs-depth budget. Use when balancing a .cool/.mcool, computing expected or P(s), making O/E matrices for compartments/loops, deciding ICE vs KR vs SCALE, choosing a resolution for a given depth, or troubleshooting NaN/all-NaN balanced matrices; route cross-sample comparison to hic-differential.
tool_type: python
primary_tool: cooler
---

## Version Compatibility

Reference examples tested with: cooler 0.10+, cooltools 0.7+, bioframe 0.7+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

cooltools standardized its API around 0.7 (functions take a `view_df` viewframe; `expected_cis` defaults to `smooth=True, aggregate_smoothed=True`). `.mcool` is multi-resolution: analysis functions take a single-resolution URI (`file.mcool::/resolutions/10000`), never the bare `.mcool`. A matrix must be balanced (a stored `weight` column) before O/E, compartments, insulation, or dots; `clr.matrix(balance=True)` on an unbalanced cooler returns all-NaN.

# Hi-C Matrix Operations

**"Make the pixels of my Hi-C matrix comparable to each other."** -> Balance (remove per-bin coverage bias under equal-visibility), then divide by distance-matched expected (remove the polymer P(s) background) to get O/E.
- Python: `cooler.balance_cooler(clr, cis_only=True, store=True)`, then `cooltools.expected_cis(clr)` and divide observed by per-diagonal expected.

## The Single Most Important Modern Insight -- Balancing Makes ONE Matrix Self-Consistent; It Does NOT Make Two Matrices Comparable

Balancing (ICE/KR) is a *within-matrix* operation: it solves for per-bin bias weights so every bin has equal genome-wide visibility, making a single map internally consistent. It does nothing to relate map A to map B. Two balanced matrices at different sequencing depth still differ in absolute magnitude, dynamic range, and noise floor -- and `rescale_marginals` makes the absolute balanced values arbitrary-scaled anyway. "I balanced both, now I'll subtract/log2-ratio them" is the single most common error in the field: the depth difference is read as biology. Cross-sample comparison requires downsampling to equal valid-pair count, distance-matched O/E, and a replicate-aware differential tool (multiHiCcompare, HiCcompare loess-over-distance, dcHiC) -- route to hic-differential.

Two corollaries that flow from the same equal-visibility model:

1. **CNV silently breaks balancing.** The premise is that every bin *should* make the same number of contacts; any deviation is technical bias. That is true for a diploid uniform-copy genome and FALSE for tumors/aneuploids -- a 3-copy region genuinely contacts ~3x more. ICE forces equal marginals and ERASES that real copy-number, then redistributes it perversely (post-ICE high-copy regions go cis-depleted, trans-enriched; Servant 2018). Use **raw counts** for CNV/SV calling (coverage IS the signal -> copy-number); plain ICE/KR on an aneuploid is a category error (use CNV-aware LOIC/CAIC for 3D structure).
2. **A balanced cis map is still dominated by distance-decay.** The A/B plaid and focal loops are a faint modulation under the P(s) background. Dividing by distance-matched expected (O/E) before eigendecomposition is mandatory, or the top eigenvector is just the decay curve, not compartments.

## Normalization-Method Taxonomy

| Method | What it does | Mechanism | When |
|--------|-------------|-----------|------|
| ICE (cooler native) | true matrix balancing | iterative proportional fitting (Sinkhorn); equalizes all marginals | default; robust, converges on sparse/low-depth where KR fails |
| KR (juicer) | true matrix balancing | Knight-Ruiz Newton solver; SAME fixed point as ICE | fast (few iterations); fails to converge on sparse/high-res maps |
| SCALE (juicer) | true matrix balancing | modern KR-family solver, more robust | juicer's default; converges where KR diverges on sparse maps |
| VC / vanilla coverage | NOT true balancing | single pass: divide by row-coverage * col-coverage (one ICE iteration) | fast robust fallback; leaves residual bias |
| VC_SQRT | NOT true balancing | divide by sqrt of coverage product (gentler than VC) | very sparse data where full balancing overfits |
| LOIC / CAIC (Servant 2018) | CNV-aware balancing | condition on copy-number; LOIC keeps the CN effect, CAIC removes it | aneuploid/tumor genomes (plain ICE is wrong here) |

KR and ICE reach the same balanced map -- choose by **convergence, not quality**: KR is faster but blows up on sparse/low-depth/high-resolution matrices; ICE is the robust default; SCALE is the juicer-side answer when KR fails.

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Need to balance a diploid map | `cooler.balance_cooler(cis_only=True, store=True)` (ICE) | robust default; cis-only is the analysis convention |
| KR failed to converge (sparse/high-res) | fall back to ICE, or SCALE on the juicer side | same fixed point; ICE/SCALE are the robust solvers |
| Tumor / aneuploid genome, 3D structure | CNV-aware LOIC/CAIC (Servant 2018) | plain ICE erases real copy-number |
| CNV / SV calling from Hi-C | RAW counts (no balancing) | coverage is the signal -> copy-number |
| Compartments at 100kb-1Mb | balance -> `expected_cis` -> O/E -> Pearson -> eigenvector | O/E removes P(s) so the plaid is visible -> compartment-analysis |
| Focal loops / TADs | balance -> `expected_cis` -> O/E | local enrichment needs the distance background removed -> loop-calling, tad-detection |
| P(s) / polymer-state diagnostic | `expected_cis(smooth=True)` -> log-derivative | the derivative reads out loop-extrusion machinery |
| Imported a juicer KR/VC weight column | check `divisive_weights` before applying | cooler weights are multiplicative, juicer's are divisive |
| Compare two conditions | downsample to equal depth, then -> hic-differential | balancing is within-matrix, not a cross-sample normalizer |
| Bare `.mcool` passed and KeyError | use `file.mcool::/resolutions/<bp>` URI | `.mcool` is a container of resolutions |

## Balance a Matrix (ICE)

**Goal:** Remove one-dimensional per-bin coverage bias so every bin has equal genome-wide visibility within this single map.

**Approach:** Mask low-coverage and blacklisted bins FIRST (mad_max on log-marginals + explicit blacklist of centromere/rDNA/unmappable), drop the first two diagonals (ligation chemistry, not 3D contact), then run cis-only ICE; the multiplicative weight vector is stored in the `weight` column.

```python
import cooler

clr = cooler.Cooler('matrix.mcool::/resolutions/10000')
bias, stats = cooler.balance_cooler(clr, cis_only=True, mad_max=5, ignore_diags=2, blacklist=None, store=True)
print('converged:', stats['converged'], 'scale:', stats['scale'])   # stats also reports var, divisive_weights

clr = cooler.Cooler('matrix.mcool::/resolutions/10000')              # re-open to see the stored weights
balanced = clr.matrix(balance=True).fetch('chr1')                    # raw[i,j] * w[i] * w[j]; masked bins -> NaN
```

`cis_only=True` is the convention for compartment/TAD/loop work -- trans signal is weak, noisy ambient ligation that pulls the bias estimates toward trans noise. `ignore_diags=2` drops the main diagonal (self-ligation/dangling ends) and first off-diagonal (undigested/religated fragments): huge untrustworthy counts that would dominate the marginals. `mad_max=5` filters bins whose **log**-marginal is >5 MAD below the median; without it a near-empty unmappable/centromeric bin gets a gigantic weight and ICE diverges. Masking (mad_max + `blacklist`) MUST precede balancing -- balancing cannot rescue a no-signal bin, it amplifies it.

CLI equivalent:

```bash
cooler balance --cis-only --mad-max 5 --ignore-diags 2 matrix.mcool::/resolutions/10000
```

## Expected: cis P(s) Curve and trans Scalar

**Goal:** Build the distance-decay background (the denominator for O/E) and the P(s) curve for diagnostics.

**Approach:** cis expected is a per-diagonal curve (one value per separation s -- this IS P(s)); trans expected is a single scalar per chromosome-pair block (trans contacts are ~distance-independent). cooltools enforces the split with two functions.

```python
import cooltools
import bioframe

clr = cooler.Cooler('matrix.mcool::/resolutions/10000')
view_df = bioframe.make_viewframe(clr.chromsizes)               # whole-chromosome regions; or arms for acrocentric genomes

cvd = cooltools.expected_cis(clr, view_df=view_df, smooth=True, aggregate_smoothed=True, ignore_diags=2)
# columns include: region1, region2, dist, dist_bp, n_valid, count.avg, balanced.avg, balanced.avg.smoothed.agg

trans_exp = cooltools.expected_trans(clr, view_df=view_df)      # one balanced.avg per region1-region2 block
```

`smooth=True` smooths P(s) in log10(distance) space (`smooth_sigma=0.1`). Pre-0.7.0 cooltools errored on raw smoothing (`clr_weight_name=None, smooth=True`, issue #456); that was fixed in 0.7.0, and raw smoothing now returns `count.avg.smoothed`. Regardless of version, balance first: a raw expected still carries per-bin coverage bias, so it is not a clean P(s)/O/E denominator.

## Observed/Expected Matrix

**Goal:** Divide out the polymer distance-decay so enrichment (loops, plaid) stands above the local background.

**Approach:** Map the per-diagonal cis expected (`balanced.avg`, keyed by `dist`) onto a dense balanced matrix by diagonal offset -- vectorized with numpy diagonal indexing, NOT an O(n^2) Python loop.

```python
import numpy as np

def oe_matrix(clr, region, cvd):
    obs = clr.matrix(balance=True).fetch(region)
    chrom = region.split(':')[0] if isinstance(region, str) else region[0]
    exp = cvd[cvd['region1'] == chrom].set_index('dist')['balanced.avg']
    exp_by_dist = exp.reindex(range(obs.shape[0])).to_numpy()              # one expected per separation s
    expected = exp_by_dist[np.abs(np.subtract.outer(np.arange(obs.shape[0]), np.arange(obs.shape[0])))]
    return obs / expected                                                  # NaN where expected is NaN (masked diags)

oe = oe_matrix(clr, 'chr1', cvd)
log_oe = np.log2(oe)                                                        # symmetric around 0 for display
```

cooltools also ships `cooltools.lib.numutils.observed_over_expected(matrix, mask)` (returns a 4-tuple `(OE, dist_bins, sum_pixels, n_pixels)`) for a self-contained dense O/E without a precomputed `cvd`.

## P(s) Log-Derivative -- the Polymer-State Diagnostic

**Goal:** Read out chromatin polymer state and loop-extrusion machinery from the shape of the contact-decay curve.

**Approach:** Take the slope of P(s) in log-log space; a reference slope near -1 over 0.1-1 Mb is the crumpled-globule background, and a loop-extrusion bump (~100kb interphase) appears as a peak in the derivative. Smooth in logspace FIRST or the derivative is pure noise.

```python
agg = cvd[(cvd['region1'] == cvd['region2']) & (cvd['dist'] > 0)].drop_duplicates('dist_bp')
slope = np.gradient(np.log(agg['balanced.avg.smoothed.agg']), np.log(agg['dist_bp']))
# slope ~ -1 over 0.1-1 Mb; a bump toward 0 near ~100kb flags cohesin loop extrusion (flattens on WAPL/RAD21 loss)
```

## Resolution-vs-Depth Budget

The achievable resolution is a function of depth and genome size, not a free choice. Rule of thumb: a bin needs ~**1000 contacts** to be reliably populated, and the number of bins scales as N^2 with the number of genomic bins -- so halving bin size quarters per-bin coverage. Coarse features are cheap, focal features are expensive:

| Feature | Resolution | Approximate depth (human) | Why |
|---------|-----------|---------------------------|-----|
| A/B compartments | 100kb-1Mb | tens of millions of valid pairs | chromosome-scale, few large bins -> cheap |
| TADs / insulation | 10-50kb | hundreds of millions | sub-Mb domains; window 5-25x the bin |
| Loops / dots | <=10kb | billions (Rao 2014 in-situ Hi-C ~ billions) | focal kb-scale pixels; coarse bins blur anchors |

Calling 10kb loops from a shallow library binned at 50kb is not a resolution choice -- there is no signal there. Choose the finest resolution where median per-bin contacts stay near ~1000.

## Per-Method Failure Modes

### Cross-sample subtraction of balanced matrices
**Trigger:** balancing two libraries then log2-ratioing/subtracting. **Mechanism:** balancing is within-matrix; balanced magnitude still scales with depth and `rescale_marginals` makes it arbitrary. **Symptom:** systematic genome-wide "differences" that track sequencing depth. **Fix:** downsample to equal valid pairs, compare O/E, use a replicate-aware tool -> hic-differential.

### Plain ICE on an aneuploid / tumor
**Trigger:** `balance_cooler` on a genome with large copy-number swings. **Mechanism:** equal-visibility forces equal marginals, erasing the real ~CN-fold coverage. **Symptom:** high-copy regions look cis-depleted/trans-enriched (Servant 2018); CNV vanishes. **Fix:** raw counts for CNV calling; LOIC/CAIC for 3D structure.

### Masking after (not before) balancing
**Trigger:** low-coverage centromere/rDNA/unmappable bins left in before ICE. **Mechanism:** a near-empty bin gets a gigantic bias weight. **Symptom:** ICE fails to converge, or stripe artifacts radiate from a few bins. **Fix:** set `mad_max` and pass `blacklist`; masking precedes balancing.

### Eigendecomposition on a balanced (non-O/E) map
**Trigger:** compartment calling skips the expected/O/E step. **Mechanism:** the distance-decay dominates the balanced cis map. **Symptom:** top eigenvector is the P(s) curve, not the A/B plaid. **Fix:** divide by `expected_cis` (O/E) before correlating/eigendecomposing.

### Crossing multiplicative and divisive weights
**Trigger:** applying an imported juicer KR/VC vector as if it were a cooler weight. **Mechanism:** cooler weights are multiplicative (raw*w*w), juicer's are divisive (raw/w/w). **Symptom:** correction inverts -- high-bias bins get MORE extreme; nothing errors. **Fix:** check the weight column's `divisive_weights` attribute before applying.

### O/E from a raw (unbalanced) expected
**Trigger:** computing P(s)/O/E from `expected_cis(clr_weight_name=None)`. **Mechanism:** raw expected still carries per-bin coverage bias, so dividing by it does not cleanly remove the polymer background; pre-0.7.0 cooltools additionally errored on `smooth=True` with raw (issue #456, fixed in 0.7.0). **Symptom:** O/E still shows coverage stripes; on old cooltools, an error demanding balanced data. **Fix:** balance first, then `expected_cis` on the balanced `weight` column.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| `ignore_diags=2` | cooler default; ICE convention | drops self-ligation/dangling (diag 0) + undigested/religated (diag 1); ligation chemistry, not 3D contact |
| `mad_max=5` | cooler default | drops bins >5 MAD below median LOG-marginal; near-empty bins otherwise get exploding weights |
| `min_nnz=10` | cooler default | <10 nonzero pixels per row is too sparse to estimate a bias reliably |
| `tol=1e-5`, `max_iters=200` | cooler defaults | convergence = variance of balanced marginals < tol; non-convergence usually = a masking problem, not a tol problem |
| `smooth_sigma=0.1` | cooltools default | Gaussian std in log10(distance) units for P(s) smoothing |
| ~1000 contacts/bin | depth-budget convention | per-bin coverage floor for a reliably populated bin; bins scale ~N^2 |
| Compartment res 100kb-1Mb | compartment scale | finer bins mix in TAD/loop structure |
| TAD/insulation res 10-50kb | domain scale | sub-Mb domains; window 5-25x the bin |
| Loop res <=10kb (needs ~billions of pairs) | Rao 2014 in-situ Hi-C | focal kb-scale; coarse bins blur loop anchors |
| P(s) slope ~ -1 over 0.1-1 Mb | crumpled/fractal globule | reference background; deviations/derivative read out polymer state |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `clr.matrix(balance=True)` all NaN | cooler not balanced | run `cooler.balance_cooler(..., store=True)` first |
| Empty / wrong-resolution result on .mcool | bare `.mcool` passed | use `file.mcool::/resolutions/<bp>` URI |
| ICE not converging | low-coverage bins not masked | raise/set `mad_max`, pass `blacklist`; do not just raise `max_iters` |
| `expected_cis(smooth=True)` errors on raw (pre-0.7.0 only) | `clr_weight_name=None` on old cooltools (issue #456, fixed 0.7.0) | upgrade to cooltools 0.7+, or balance first / `smooth=False` |
| O/E enrichment inverted on a tumor | ICE applied to an aneuploid | use raw counts / LOIC/CAIC; equal-visibility is violated |
| Imported weight makes bias worse | divisive juicer weight applied as multiplicative | check `divisive_weights`; reciprocate if needed |
| Cross-condition "difference" tracks depth | subtracting balanced matrices | downsample + O/E + replicate-aware test -> hic-differential |

## References

- Imakaev et al. 2012 *Nat Methods* 9:999-1003 -- ICE iterative correction.
- Knight & Ruiz 2013 *IMA J Numer Anal* 33(3):1029-1047 -- KR matrix-balancing algorithm.
- Rao et al. 2014 *Cell* 159(7):1665-1680 -- in-situ Hi-C, KR norm, VC/VC_SQRT, kilobase loops (depth budget).
- Cournac et al. 2012 *BMC Genomics* 13:436 -- sequential/vanilla-coverage (SCN/VC) normalization.
- Servant et al. 2018 *BMC Bioinformatics* 19:313 -- CNV-aware normalization (LOIC/CAIC); ICE erases copy-number.
- Abdennur & Mirny 2020 *Bioinformatics* 36(1):311-316 -- cooler.
- Open2C, Abdennur et al. 2024 *PLoS Comput Biol* 20(5):e1012067 -- cooltools.
- Open2C, Abdennur et al. 2024 *Bioinformatics* 40(2):btae088 -- bioframe.

## Related Skills

- hic-data-io - Load the cooler files this skill balances; divisive-vs-multiplicative weight naming
- compartment-analysis - Consumes the O/E this skill produces for eigenvector calling
- tad-detection - Insulation needs a cis-balanced matrix
- loop-calling - Dots need balanced + expected as prerequisites
- hic-differential - Cross-sample comparison; the right home for subtracting/ratioing conditions
- hic-visualization - Render balanced/O/E/log matrices
- copy-number/cnv-visualization - Raw-count CNV from Hi-C when balancing would erase copy-number
- genome-intervals/bigwig-tracks - Export the P(s)/expected or eigenvector as a bigWig

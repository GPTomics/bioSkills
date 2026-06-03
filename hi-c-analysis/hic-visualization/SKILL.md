---
name: bio-hi-c-analysis-hic-visualization
description: Renders Hi-C contact matrices honestly and reproducibly with matplotlib, cooltools, HiCExplorer, pyGenomeTracks, FAN-C, CoolBox, and plotgardener. Covers the raw/ICE-balanced/observed-over-expected transform choice, LogNorm vs symmetric-diverging colormaps with vmax/percentile clipping, resolution-to-feature matching (compartments 100-500kb, TADs 10-40kb, loops 5-10kb), square vs rotated-triangle track-stacking, NaN/white-stripe handling, virtual 4C, APA/saddle/on-diagonal pileups, two-condition side-by-side and log2-ratio maps, and interactive (HiGlass) vs scripted-static publication figures. Use when plotting a contact matrix, choosing a normalization or color scale, building a multi-track Hi-C figure, making a virtual 4C profile, piling up loops/boundaries, or comparing two conditions.
tool_type: python
primary_tool: matplotlib
---

## Version Compatibility

Reference examples tested with: cooler 0.10+, cooltools 0.7+, matplotlib 3.8+, bioframe 0.7+, HiCExplorer 3.7+, pyGenomeTracks 3.9+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

Notes specific to this skill: `.mcool` is multi-resolution -- pass a single-resolution URI (`file.mcool::/resolutions/10000`), never the bare path. A cooler must be balanced (`cooler balance`) before `matrix(balance=True)` returns anything but NaN. `cooltools.pileup` returns a stack of shape `(n_features, D, D)` -- aggregate over `axis=0`. cooltools standardised on `view_df`/`expected_df` arguments around 0.7+; verify with `help(cooltools.pileup)` before chaining.

# Hi-C Visualization

**"Plot my Hi-C contact matrix"** -> Choose a transform (raw / ICE-balanced / observed-over-expected), a matched colormap+norm (LogNorm for counts, symmetric-diverging for O/E and ratios), and a resolution that fits the feature; render as a square map or a rotated triangle for track-stacking, with NaN bins shown explicitly.
- Python: `clr.matrix(balance=True).fetch(region)` then `ax.matshow(m, norm=LogNorm(...), cmap='fall')`
- CLI: `hicPlotMatrix --matrix m.cool --region chr1:50-60Mb --log1p --colorMap fall -o out.png`

## The Single Most Important Modern Insight -- The Colorscale Is Where Hi-C Figures Lie

The same matrix under raw / ICE-balanced / observed-over-expected tells three *different biological stories*, and the choice is not cosmetic -- it decides which biology is legible. A balanced map is still dominated by the polymer distance-decay gradient (the bright diagonal falling off as ~P(s)); it shows TADs but washes out compartments and loops. Dividing by the distance-matched expected and taking `log2(O/E)` with a SYMMETRIC diverging cmap (`coolwarm`/`RdBu_r`, `vmin=-vmax`) removes that gradient and makes the compartment checkerboard and loop corner-dots suddenly visible -- they were always in the data. The corollary is a reviewer's reflex: a "no compartments / no loops" claim plotted on a balanced (not O/E) map is unsupportable. A reviewer-grade figure is one that can be reconstructed from the legend -- it states (1) the normalization (raw / ICE-balanced / O/E / log2-ratio), (2) the color scale (LogNorm vs symmetric-diverging) with its limits or clip percentile, and (3) the bin resolution. If those three are absent, the figure is neither interpretable nor reproducible.

## Transform Taxonomy

| Transform | Norm + cmap | What it shows | When |
|-----------|-------------|---------------|------|
| Raw counts | `LogNorm`, sequential (`fall`) | depth + per-bin coverage bias; white stripes are artifacts | QC sanity check only -- almost never the science figure |
| ICE-balanced | `LogNorm` vmin~1e-4..1e-1, `fall` | TADs + the distance-decay gradient; loops/compartments washed out | the honest "raw structure" map; track-stacking context |
| Observed/Expected | `log2`, symmetric `coolwarm`/`RdBu_r`, `vmin=-vmax` | compartment checkerboard + loop corner-dots | compartments, loops, any focal-enrichment claim |
| log2(cond1/cond2) | symmetric `RdBu_r`, `vmin=-vmax`, white=0 | gained/lost contacts | two conditions (balance + depth-match FIRST) |

## Layout Taxonomy

| Layout | Tool | Mechanism | When |
|--------|------|-----------|------|
| Square map | matplotlib `matshow`/`pcolormesh`, FAN-C `HicPlot2D` | symmetric 2D matrix | matrix itself is the result; inter-region rectangle; difference map |
| Rotated triangle | pyGenomeTracks/HiCExplorer `hic_matrix`, FAN-C `HicPlot`, plotgardener `plotHicTriangle`, CoolBox `style='triangular'` | 45deg shear, keep upper half, diagonal on top | STACKING genome-browser tracks below on a shared x-axis |
| Pileup (APA/saddle/on-diagonal) | `cooltools.pileup`/`saddle`, coolpup.py | average snippets over a feature set | the only honest genome-wide claim from sparse data |
| Virtual 4C | one matrix row, FAN-C `HicSlicePlot` | 1D profile from a viewpoint bin | compare against a real 4C anchor |
| Interactive | HiGlass | multires `.mcool` pan/zoom | exploration -- NOT a reproducible figure |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Show compartments / loops | `log2(O/E)`, symmetric `coolwarm`, `vmin=-vmax` | balanced map's gradient hides them |
| Show TADs / domains | balanced `LogNorm`, `fall`, 10-40kb | domain insulation lives at sub-Mb scale |
| Matrix + genes + ChIP + insulation stack | -> data-visualization/genome-tracks (pyGenomeTracks `hic_matrix`) | config = reproducible provenance; library does the shear |
| Quantify compartment strength | saddle plot (phase E1 first) -> compartment-analysis | corners give the single strength number |
| Validate a loop SET genome-wide | APA pileup, `log2(O/E)`, symmetric | one loop is invisible; 10k averaged is solid |
| Compare two conditions | side-by-side same-scale OR log2-ratio | balance + depth-match both FIRST |
| Export eigenvector / insulation as a track | -> genome-intervals/bigwig-tracks | bigWig feeds the track stack |
| Explore to find a region/resolution | HiGlass, then reproduce in a scripted tool | interactive != publication |

## Square Contact Map (Balanced, Log-Scaled, NaN Shown)

**Goal:** Render a balanced cis matrix for one region with honest dynamic range and visible masked bins.

**Approach:** Fetch the balanced matrix at a single-resolution URI, set `vmax` from a high off-diagonal percentile (report it), use `LogNorm`, and explicitly color NaN bins with `set_bad` so masked regions read as gray rather than as "zero contact".

```python
import cooler, numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import cooltools.lib.plotting   # registers the 'fall' cmap; needs matplotlib < 3.9 with cooltools 0.7.x (else use a stock cmap like 'afmhot_r')

clr = cooler.Cooler('matrix.mcool::/resolutions/10000')
region = ('chr1', 50_000_000, 60_000_000)
m = clr.matrix(balance=True).fetch(region)

vmax = np.nanpercentile(m[m > 0], 99.5)   # report this percentile in the legend
cmap = plt.get_cmap('fall').copy(); cmap.set_bad('lightgray')   # NaN bins shown, not white
fig, ax = plt.subplots(figsize=(7, 7))
im = ax.matshow(m, norm=LogNorm(vmin=vmax * 1e-3, vmax=vmax), cmap=cmap)
fig.colorbar(im, ax=ax, fraction=0.046, label='balanced (ICE)')
```

## Observed/Expected Divergent Map

**Goal:** Make the compartment checkerboard and loop corner-dots legible by removing the polymer distance-decay background.

**Approach:** Compute the cis expected with cooltools, fetch the matched-distance expected per pixel, divide observed by expected and take `log2`, then plot with a symmetric diverging cmap centered at 0 -- asymmetric limits move the white midpoint off zero and make neutral regions read as enriched.

```python
import cooltools, bioframe

view_df = bioframe.make_viewframe(clr.chromsizes)
expected = cooltools.expected_cis(clr, view_df=view_df, nproc=4)

chrom = region[0]
exp_by_diag = expected.query('region1 == @chrom')['balanced.avg'].to_numpy()   # expected per genomic separation
m = clr.matrix(balance=True).fetch(region)
i, j = np.indices(m.shape)
oe_mtx = m / exp_by_diag[np.abs(i - j)]   # divide each pixel by its distance-matched expected

v = 2.0   # symmetric clip; |log2(O/E)| up to ~2 is the usual readable range
fig, ax = plt.subplots(figsize=(7, 7))
im = ax.matshow(np.log2(oe_mtx), cmap='coolwarm', vmin=-v, vmax=v)   # vmin=-vmax mandatory
fig.colorbar(im, ax=ax, fraction=0.046, label='log2(obs/exp)')
```

`expected['balanced.avg']` is the per-diagonal expected; indexing it by `|i-j|` broadcasts it to a full per-pixel expected matrix.

## Rotated Triangle for Track-Stacking

**Goal:** Hang the contact map above aligned genome-browser tracks (genes, ChIP, insulation) on a shared x-axis.

**Approach:** Prefer a library that owns the 45deg shear and the `depth` crop -- pyGenomeTracks/HiCExplorer (`file_type = hic_matrix`), FAN-C, plotgardener, or CoolBox -- because hand-rolling the `Affine2D` shear is where `extent`/`aspect` alignment bugs live. The matplotlib reference below is for a single panel; for a real stack, route to data-visualization/genome-tracks.

```bash
# HiCExplorer / pyGenomeTracks: config-driven, reproducible. depth = how far up the diagonal.
# A 2 Mb TAD needs depth >= ~2_000_000 or it is silently truncated.
hicPlotTADs --tracks tracks.ini --region chr1:50000000-60000000 -o stack.png
```

```python
from matplotlib.transforms import Affine2D
# matplotlib single-panel reference: shear the square map onto the diagonal.
t = Affine2D().rotate_deg(45) + ax.transData
im = ax.pcolormesh(np.log2(oe_mtx), cmap='coolwarm', vmin=-v, vmax=v)
im.set_transform(t)
ax.set_ylim(0, m.shape[0])   # crop the depth; the y-axis is genomic SEPARATION, not a 2nd coordinate
```

## Virtual 4C from a Viewpoint

**Goal:** Extract a 1D contact profile from one viewpoint bin to compare against a real 4C experiment.

**Approach:** Take the viewpoint row from the balanced chromosome matrix; the near-cis distance-decay spike swamps distal signal on a linear axis, so plot on log-y (or mask the +/- few bins around the viewpoint) and say which. Cross-condition profiles must be balanced, depth-matched, and on identical y-axes.

```python
res = clr.binsize
vp_bin = (55_000_000 // res) - (50_000_000 // res)   # viewpoint index within the region
profile = clr.matrix(balance=True).fetch(region)[vp_bin, :]
fig, ax = plt.subplots(figsize=(11, 2.5))
ax.semilogy(np.arange(len(profile)) * res / 1e6 + 50, profile)   # log-y: the near-cis spike lies on linear
ax.axvline(55, color='red', ls='--')
```

## APA Pileup over a Loop Set

**Goal:** Validate a loop call set genome-wide by averaging the contact signal centered on every loop's anchor pair.

**Approach:** Pass BEDPE features and the cis expected to `cooltools.pileup` for an observed/expected stack, average over the feature axis (`axis=0`), and plot `log2` with a symmetric cmap; the center pixel is the loop, the APA score is center / a corner-background patch (Rao 2014 lower-left 3x3 convention).

```python
import pandas as pd
loops = pd.read_csv('loops.bedpe', sep='\t')   # chrom1,start1,end1,chrom2,start2,end2
stack = cooltools.pileup(clr, loops, view_df=view_df, expected_df=expected, flank=100_000)
apa = np.nanmean(stack, axis=0)   # stack is (n_features, D, D) -> average over features
c = apa.shape[0] // 2
apa_score = apa[c, c] / np.nanmean(apa[-3:, :3])   # center / lower-left 3x3 background
fig, ax = plt.subplots(figsize=(5, 5))
im = ax.matshow(np.log2(apa), cmap='coolwarm', vmin=-1, vmax=1)
ax.set_title(f'APA score {apa_score:.2f}')
```

For on-diagonal pileups over CTCF sites, strand-orient before averaging (`stack[mask] = stack[mask][:, ::-1, ::-1]` for `-` strand) or convergent/divergent signals cancel. For HiChIP/PLAC-seq anchored loops, route to loop-calling and the peak-anchored pileup conventions there.

## Two-Condition Comparison

**Goal:** Show a contact change between two conditions without it being a depth/coverage artifact.

**Approach:** Both matrices must be ICE-balanced AND depth-matched (downsample the deeper library to equal valid pairs) BEFORE ratioing. Then either side-by-side panels on an IDENTICAL cmap/norm/vmin/vmax/resolution, or a single `log2(cond1/cond2)` divergent map with symmetric limits and white = no change; grey out very-distal noise-amplified bins.

```python
m1 = clr1.matrix(balance=True).fetch(region)
m2 = clr2.matrix(balance=True).fetch(region)   # clr1, clr2 already depth-equalized upstream
ratio = np.log2((m1 + 1e-5) / (m2 + 1e-5))   # pseudocount tames divide-by-small off-diagonal
fig, ax = plt.subplots(figsize=(7, 7))
im = ax.matshow(ratio, cmap='RdBu_r', vmin=-2, vmax=2)   # symmetric, white=no change
fig.colorbar(im, ax=ax, fraction=0.046, label='log2(cond1/cond2)')
```

Quantitative replicate-aware differential testing (not just a figure) lives in hic-differential.

## Per-Method Failure Modes

### Negative claim on a balanced map
**Trigger:** "no compartments/loops" read off a balanced (not O/E) map. **Mechanism:** the polymer distance-decay gradient dominates balanced data and hides checkerboard/dots. **Symptom:** features absent that O/E would reveal. **Fix:** replot as `log2(O/E)` with a symmetric cmap before making any negative claim.

### Asymmetric divergent limits
**Trigger:** `vmin != -vmax` on an O/E or log2-ratio map. **Mechanism:** zero (no change/enrichment) is no longer the white midpoint. **Symptom:** neutral regions read as enriched or depleted; reviewers flag it. **Fix:** `vmin=-v, vmax=v` (or `TwoSlopeNorm(vcenter=0)`).

### vmax = data max
**Trigger:** letting the heavy-tailed diagonal set vmax. **Mechanism:** a few super-bins are orders of magnitude above the bulk. **Symptom:** the whole map looks empty/dark; over-clipping the other way fabricates structure. **Fix:** vmax = a stated high off-diagonal percentile (95th-99.5th).

### Resolution mismatched to feature
**Trigger:** plotting at whatever the `.mcool` defaults to. **Mechanism:** loops at 100kb are averaged away (oversmoothing); compartments at 5-10kb mix in TAD/loop noise. **Symptom:** vanished dots or a noisy checkerboard. **Fix:** compartments 100-500kb, TADs 10-40kb, loops 5-10kb (Micro-C 1-2kb).

### Interpolated matrix quantified
**Trigger:** `interp_nan`/`adaptive_coarsegrain` (or scHi-C smoothing) then measuring on the result. **Mechanism:** interpolation/imputation fabricates contacts for DISPLAY. **Symptom:** a filled centromere looks like contiguous chromatin; "structure" that is the smoother's prior. **Fix:** fill for display only; quantify on the raw/balanced matrix and disclose the smoother.

### Triangle depth too small
**Trigger:** `depth` < the largest feature. **Mechanism:** the triangle crop truncates separations above `depth`. **Symptom:** a 2 Mb TAD silently cut off. **Fix:** set `depth >= ~feature size`; remember the y-axis is genomic separation.

### Unequal-depth ratio
**Trigger:** `log2(cond1/cond2)` on un-depth-matched or unbalanced maps. **Mechanism:** a global depth difference is a uniform multiplicative offset. **Symptom:** a whole-map color shift read as biology. **Fix:** ICE-balance and downsample to equal valid pairs first.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Compartment resolution 100-500kb | compartment scale (Lieberman-Aiden 2009) | checkerboard is Mb-scale; finer bins add noise, not detail |
| TAD resolution 10-40kb | domain scale (Dixon 2012) | insulation/boundary structure lives at sub-Mb |
| Loop resolution 5-10kb (Micro-C 1-2kb) | focal-contact scale (Rao 2014) | a loop is a ~10kb focal pixel; coarse bins blur it, too-fine buries it in Poisson noise |
| vmax = 95th-99.5th off-diagonal percentile | heavy-tailed counts | data-max vmax leaves the map dark; report the percentile |
| Divergent limits symmetric `vmin=-vmax` | zero must be the midpoint | asymmetric limits misplace the white neutral point |
| APA flank +/- 100kb | corner-background convention | too small contaminates the corner; too large averages in neighbors |
| APA score = center / lower-left 3x3 | Rao 2014 | standard center-to-background loop enrichment ratio |
| HiGlass zoom levels < ~5x apart | Kerpedjiev 2018 | adjacent resolutions must be close for smooth multires rendering |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `matrix(balance=True)` all NaN | cooler not balanced | run `cooler balance` / `cooler.balance_cooler` first |
| Empty / wrong-resolution result | bare `.mcool` passed | use `file.mcool::/resolutions/<bp>` |
| White stripes mistaken for "no contact" | NaN bins left at default | `cmap.set_bad('lightgray')` to render masked bins |
| `cmap='fall'` KeyError | colormap not registered | `import cooltools.lib.plotting` first |
| `ImportError` on `import cooltools.lib.plotting` | matplotlib >= 3.9 dropped `register_cmap` (cooltools 0.7.x) | pin matplotlib < 3.9, or use a stock cmap (`'afmhot_r'`) |
| Pileup looks averaged-out / scrambled | wrong nanmean axis | aggregate over `axis=0` (stack is `(n_features, D, D)`) |
| Empty region / no overlap | chrom naming (`chr1` vs `1`) | harmonize names across cooler, BED/BEDPE, fasta |
| `AttributeError` on cooltools fn | pre-0.7 vs 0.7+ API | `help(cooltools.<fn>)`; update to the `view_df`/`expected_df` signature |

## References

- cooler: Abdennur N, Mirny LA. Cooler: scalable storage for Hi-C data and other genomically labeled arrays. *Bioinformatics* 2020;36(1):311-316.
- cooltools: Open2C, Abdennur N, Abraham S, Fudenberg G, et al. Cooltools: enabling high-resolution Hi-C analysis in Python. *PLoS Comput Biol* 2024;20(5):e1012067.
- HiCExplorer: Ramirez F, Bhardwaj V, Arrigoni L, et al. High-resolution TADs reveal DNA sequences underlying genome organization in flies. *Nat Commun* 2018;9(1):189.
- pyGenomeTracks: Lopez-Delisle L, Rabbani L, Wolff J, et al. pyGenomeTracks: reproducible plots for multivariate genomic datasets. *Bioinformatics* 2021;37(3):422-423.
- FAN-C: Kruse K, Hug CB, Vaquerizas JM. FAN-C: a feature-rich framework for the analysis and visualisation of chromosome conformation capture data. *Genome Biol* 2020;21(1):303.
- CoolBox: Xu W, Zhong Q, Lin D, et al. CoolBox: a flexible toolkit for visual analysis of genomics data. *BMC Bioinformatics* 2021;22(1):489.
- plotgardener: Kramer NE, Davis ES, Wenger CD, et al. Plotgardener: cultivating precise multi-panel figures in R. *Bioinformatics* 2022;38(7):2042-2045.
- HiGlass: Kerpedjiev P, Abdennur N, Lekschas F, et al. HiGlass: web-based visual exploration and analysis of genome interaction maps. *Genome Biol* 2018;19(1):125.
- coolpup.py: Flyamer IM, Illingworth RS, Bickmore WA. Coolpup.py: versatile pile-up analysis of Hi-C data. *Bioinformatics* 2020;36(10):2980-2985.
- A/B compartments: Lieberman-Aiden E, van Berkum NL, Williams L, et al. Comprehensive mapping of long-range interactions reveals folding principles of the human genome. *Science* 2009;326(5950):289-293.
- TAD directionality index: Dixon JR, Selvaraj S, Yue F, et al. Topological domains in mammalian genomes identified by analysis of chromatin interactions. *Nature* 2012;485(7398):376-380.
- Loops/HiCCUPS/APA: Rao SSP, Huntley MH, Durand NC, et al. A 3D map of the human genome at kilobase resolution reveals principles of chromatin looping. *Cell* 2014;159(7):1665-1680.

## Related Skills

- hic-data-io - Load the cooler/.mcool files and zoomify for multires HiGlass tilesets
- matrix-operations - Balancing and O/E that the divergent map depends on
- compartment-analysis - Eigenvector phasing behind the saddle plot
- tad-detection - Insulation/boundary tracks stacked under the triangle
- loop-calling - Loop calls and peak-anchored pileup conventions visualized here
- hic-differential - Replicate-aware testing behind the two-condition comparison
- data-visualization/genome-tracks - Config-driven multi-track stacks (pyGenomeTracks hic_matrix)
- genome-intervals/bigwig-tracks - Export eigenvector/insulation as bigWig for the track stack
- data-visualization/heatmaps-clustering - General heatmap color/norm conventions

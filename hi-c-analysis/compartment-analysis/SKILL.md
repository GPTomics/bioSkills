---
name: bio-hi-c-analysis-compartment-analysis
description: Detects A/B chromatin compartments from balanced Hi-C contact matrices via eigenvector decomposition of the distance-normalized, Pearson-correlated cis matrix with cooltools (eigs_cis), then orients (phases) the compartment eigenvector against a GC or gene-density track so the active (A) sign is not arbitrary. Covers the eigenvector-is-a-choice problem (per-arm view_df to remove the centromere gradient; picking the eigenvector by max correlation with activity, not by eigenvalue), GC phasing with bioframe.frac_gc, resolution choice (100kb-1Mb), saddle plots and saddle_strength for compartmentalization strength, the cohesin-loss-strengthens-compartments result, subcompartments (SNIPER/Calder/dcHiC), and cross-condition compartment switching. Use when calling A/B compartments, computing E1/eigenvectors, phasing the eigenvector, building saddle plots, choosing a compartment resolution, quantifying compartment strength, or comparing compartmentalization across conditions.
tool_type: python
primary_tool: cooltools
---

## Version Compatibility

Reference examples tested with: cooler 0.10+, cooltools 0.7+, bioframe 0.7+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

cooltools had a major API shift around 0.5 -> 0.7+ (functions standardized on `view_df`/viewframe arguments; `eigs_cis`, `expected_cis`, `saddle` signatures changed). The cooler MUST be balanced before any compartment analysis: `clr.matrix(balance=True)` requires a stored `weight` column. A `.mcool` is multi-resolution -- pass a single-resolution URI (`file.mcool::/resolutions/100000`), not the bare `.mcool`. The `phasing_track` MUST share the cooler's exact binning or phasing silently no-ops. If code throws ImportError, AttributeError, or TypeError, introspect the installed package and adapt the example to match the actual API rather than retrying.

# A/B Compartment Analysis

**"Which regions of my genome are in the active vs inactive compartment?"** -> Distance-normalize the cis matrix, take an eigenvector of its Pearson correlation matrix, then orient it by GC/gene density so positive = A (active) -- but verify the kept eigenvector is the compartment one, not an arm gradient.
- Python: `cooltools.eigs_cis(clr, gc_track, view_df=arms, n_eigs=3, sort_metric='pearsonr')`

## The Single Most Important Modern Insight -- E1 Is a Choice, Not an Output, and Its Sign Is Arbitrary Until Phased

The two most damaging beginner assumptions are "E1 = compartments" and "positive E1 = active." Both are false out of the box, and both fail *silently* -- the pipeline runs, returns a track, and is wrong.

1. **E1 is not guaranteed to be the compartment track.** On a *whole-chromosome* O/E correlation matrix the largest eigenvalue very often belongs to a smooth p-arm-vs-q-arm or centromere-to-telomere GRADIENT, not the plaid A/B checkerboard -- the real compartment signal then lands in E2 or E3. cooltools' own docs concede the first eigenvector "occasionally describes chromosomal arms or translocation blowouts." The compartment eigenvector is the one with the largest |correlation| to an activity track (GC, gene density, H3K27ac), not the one with the largest eigenvalue. The structural fix removes the gradient at the source: run `eigs_cis` per chromosome ARM (a `view_df` split at centromeres, from `bioframe.make_chromarms`), so the arm gradient is never in the within-arm matrix. Set `sort_metric='pearsonr'` so the returned eigenvectors are ordered by GC correlation, not eigenvalue -- otherwise the arm gradient is reported as "E1." A monotonic "compartment track" with no sign flips across a chromosome is the failure signature of a captured arm gradient.

2. **The sign is arbitrary until phased.** Eigenvectors are defined up to sign; the positive lobe is meaningless and can differ per chromosome AND per sample. The eigenvector MUST be oriented with an external active-chromatin track via the `phasing_track` argument so A = positive. GC content is the field default (it needs no extra assay and tracks compartment A; Lieberman-Aiden 2009 *Science* 326:289) -- compute it with `bioframe.frac_gc` at the compartment resolution, exactly matching the cooler's binning. Wrong/weak phasing flips A<->B silently, and every downstream saddle, switch call, and differential result inverts with no error. This is a classic source of irreproducible compartment papers.

3. **Compartments are an equilibrium phenomenon decoupled from TADs/loops.** Compartments = microphase separation of A/B chromatin states (cohesin-independent; survive cohesin loss, Schwarzer 2017 *Nature* 551:51, and CTCF loss, Nora 2017 *Cell* 169:930). TADs/loops = ATP-driven loop extrusion stalled at CTCF. Removing cohesin reinforces compartments while erasing TADs (Schwarzer 2017 reports reinforced compartmentalization on Nipbl loss; Rao 2017 *Cell* 171:305 eliminates all loop domains with compartments retained) -- loop extrusion actively mixes chromatin across compartment boundaries, so removing the extruder lets microphase separation run to completion (Nuebler 2018 *PNAS* 115:E6697). A preserved-or-stronger saddle after a cohesin/Nipbl/RAD21 perturbation is the EXPECTED result, not a bug; compartment-strength and TAD-strength are antagonistic. If a CTCF/cohesin perturbation makes compartments *vanish*, suspect a phasing artifact, not biology.

## Method / Output Taxonomy

| Output | Tool / call | What it is | When |
|--------|-------------|-----------|------|
| A/B eigenvector (E1) | `cooltools.eigs_cis` (cis, per-arm) | leading GC-phased eigenvector of the cis O/E correlation matrix; sign = A/B | standard A/B call, single map, per chromosome arm |
| Genome-wide A/B | `cooltools.eigs_trans` | eigenvector of inter-chromosomal blocks; immune to the cis arm-gradient | whole-genome A/B consensus with deep trans coverage |
| Compartment strength | `cooltools.saddle` + `saddle_strength` | (AA+BB)/(AB+BA) corner ratio of the saddle | comparing compartmentalization across conditions |
| 5-6 subcompartments | SNIPER (Xiong & Ma 2019 *Nat Commun* 10:5069) | autoencoder imputes inter-chr contacts -> MLP classifies A1/A2/B1/B2/B3 at 100kb | deep inter-chr data; Rao-style subcompartments |
| Continuous compartment rank | Calder (Liu 2021 *Nat Commun* 12:2439) | intra-chr divisive hierarchical clustering -> 0-1 multi-scale rank | cross-cell-line repositioning; modest coverage |
| Differential compartments | dcHiC (Chakraborty 2022 *Nat Commun* 13:6827) | quantile-normalized scores + multivariate Mahalanobis distance + significance; solves cross-sample sign flips | >=2 samples, "which bins switch A<->B" |
| Single-cell compartment | scA/B (Tan 2018 Dip-C *Science*) | CpG/activity proxy per locus; do NOT eigendecompose one sparse cell | scHi-C, ~20-50k contacts/cell |

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Matrix not yet balanced | `cooler balance` first (-> matrix-operations) | unbalanced -> O/E is all-NaN; meaningless eigenvector |
| Standard A/B call, one map | `eigs_cis` per arm at 100kb-1Mb, phase by GC | compartments are chromosome-scale; arms remove the centromere gradient |
| E1 looks monotonic / no sign flips | inspect E2/E3, pick by max |corr| with GC; or split by arm | E1 captured an arm/translocation gradient, not compartments |
| Sign of A/B seems inverted | confirm `phasing_track` is at the cooler's binning | weak/mismatched phasing flips A<->B silently |
| Want compartment STRENGTH | `saddle` + `saddle_strength`, fixed extent across samples | a single eigenvector does not quantify strength |
| Want 5-6 subcompartments | SNIPER or Calder, NOT more eigenvectors | subcompartments need inter-chr ML or hierarchical clustering, not `n_eigs` |
| Two+ conditions, compartment shift | -> hic-differential (dcHiC) | replicate-aware, sign-coherent across the cohort; hand-diffing eigenvectors flips signs |
| Single-cell Hi-C | scA/B (Dip-C), not a per-cell eigenvector | one sparse cell is too noisy to eigendecompose |
| Annotate switched bins with marks | -> chip-seq/chromatin-state-segmentation, chip-seq/peak-annotation | overlay ChromHMM/histone state on compartment calls |
| Render the eigenvector/saddle | -> hic-visualization; export bigWig -> genome-intervals/bigwig-tracks | track/heatmap conventions live there |

## Per-Arm Eigenvector with GC Phasing

**Goal:** Assign each genomic bin to the active (A) or inactive (B) compartment with a non-arbitrary sign, avoiding the centromere arm-gradient artifact.

**Approach:** Build a per-arm `view_df` (split at centromeres) so the arm gradient never enters the matrix; compute a GC-content phasing track at the cooler's exact binning; run `eigs_cis` with the GC track and `sort_metric='pearsonr'` so eigenvectors are ordered by GC correlation; then take the GC-correlated eigenvector as the compartment track.

```python
import cooler
import cooltools
import bioframe

clr = cooler.Cooler('matrix.mcool::/resolutions/100000')   # 100kb: compartments are coarse-scale
chromsizes = clr.chromsizes
cens = bioframe.fetch_centromeres('hg38')
arms = bioframe.make_chromarms(chromsizes, cens)            # per-arm view removes the centromere gradient
arms = arms[arms.chrom.isin(clr.chromnames)].reset_index(drop=True)

genome = bioframe.load_fasta('hg38.fa')                     # FASTA index (.fai) must exist
bins = clr.bins()[:][['chrom', 'start', 'end']]
gc = bioframe.frac_gc(bins, genome)                         # phasing track at the cooler's exact binning

eigvals, eigvecs = cooltools.eigs_cis(clr, gc, view_df=arms, n_eigs=3, sort_metric='pearsonr')
eigvecs['compartment'] = ['A' if e > 0 else 'B' for e in eigvecs['E1']]   # GC-phased: positive E1 = A
```

After phasing, sanity-check that `E1` correlates with `gc['GC']` (sign and magnitude). If the strongest correlation is in `E2`/`E3`, that component -- not `E1` -- is the compartment track; re-derive the call from it.

## Compartment Strength via Saddle Plot

**Goal:** Quantify how strongly the genome demixes into A and B with a single comparable number across conditions.

**Approach:** Compute the distance-decay expected; pass the cooler, the expected, and the *phased* E1 eigenvector to `saddle`, which digitizes E1 into quantile groups internally (via `qrange`) and aggregates O/E into a 2D table of same-vs-cross-compartment interactions; then read `saddle_strength` (the (AA+BB)/(AB+BA) corner ratio) at one fixed extent, applied identically to every sample compared.

```python
N_GROUPS = 38            # quantile groups for digitizing E1; ~30-50 is conventional (cooltools tutorial)
Q_LO, Q_HI = 0.025, 0.975   # trim the extreme 2.5% tails before digitizing to resist outlier bins

expected = cooltools.expected_cis(clr, view_df=arms)        # has the 'balanced.avg' column saddle needs
track = eigvecs[['chrom', 'start', 'end', 'E1']]           # the SAME phased E1 used for the A/B call
interaction_sum, interaction_count = cooltools.saddle(
    clr, expected, track, 'cis', n_bins=N_GROUPS, qrange=(Q_LO, Q_HI), view_df=arms
)
strength = cooltools.api.saddle.saddle_strength(interaction_sum, interaction_count)   # 1D array; lives in cooltools.api.saddle, not top level
EXTENT = N_GROUPS // 5   # read strength at the top/bottom ~20% of bins; pick one extent, use it everywhere
score = strength[EXTENT]
```

`saddle_strength` returns an ARRAY (cumulative corner ratio over increasing extent), not a scalar -- there is no canonical single number, so choose an extent and apply it identically across compared samples. Mismatched `n_bins`, resolution, `qrange`, or extent make strengths incomparable. Remember: a preserved-or-*higher* strength after cohesin/Nipbl/RAD21 loss is the expected result.

## Compare Compartments Across Conditions

**Goal:** Find bins that switch A<->B between two conditions without being fooled by per-sample sign flips.

**Approach:** Do NOT independently phase two eigenvectors and diff them bin-by-bin -- a weak per-chromosome GC correlation can flip the sign in one sample only, manufacturing fake "switches." Use dcHiC, which computes eigenvectors on quantile-normalized scores in a shared framework (sign-coherent across the cohort) and reports a multivariate significance per bin. Route this to hic-differential.

```python
# Hand-diffing is only safe when you have CONFIRMED both eigenvectors are sign-coherent (same arms phased
# to the same GC track with strong correlation). Otherwise use dcHiC -- see hic-differential.
import pandas as pd
merged = eig1.merge(eig2, on=['chrom', 'start', 'end'], suffixes=('_1', '_2'))
merged['switch'] = (merged['E1_1'] > 0) != (merged['E1_2'] > 0)   # only meaningful if both are phased coherently
```

## Per-Method Failure Modes

### Eigenvector captured the arm gradient
**Trigger:** `eigs_cis` run per whole chromosome (no per-arm `view_df`) and/or `sort_metric=None`. **Mechanism:** the largest eigenvalue belongs to the smooth p-vs-q arm / centromere gradient, not the A/B checkerboard. **Symptom:** a monotonic "compartment track" across a chromosome with no sign flips; weak correlation of E1 with GC. **Fix:** run per chromosome arm (`bioframe.make_chromarms`); set `sort_metric='pearsonr'`; pick the eigenvector with the largest |corr| to GC.

### Eigenvector sign not phased
**Trigger:** `eigs_cis` called with `phasing_track=None`. **Mechanism:** the sign of an eigenvector is mathematically arbitrary. **Symptom:** active euchromatin lands in "B"; A/B inverted relative to GC; per-chromosome sign inconsistency. **Fix:** pass a GC (or gene-density / H3K27ac) `phasing_track` so positive E1 = A.

### Phasing track at the wrong binning
**Trigger:** GC/activity track computed at a different resolution than the cooler. **Mechanism:** cooltools aligns the track to the cooler bins; a mismatch yields garbage correlations or a silent no-op. **Symptom:** phasing has no effect, or signs are random. **Fix:** compute the track on `clr.bins()` at the exact compartment resolution.

### Calling compartments at TAD/loop resolution
**Trigger:** `eigs_cis` at 5-25kb. **Mechanism:** compartments are a 100kb-1Mb feature; fine bins are sparse and dominated by TAD/loop structure and noise. **Symptom:** a noisy, jagged E1 that does not correlate with GC. **Fix:** call at 100kb-1Mb (250kb common; up to 1Mb for shallow data).

### Expecting subcompartments from more eigenvectors
**Trigger:** raising `n_eigs` to "get A1/A2/B1/B2/B3". **Mechanism:** Rao's 6 subcompartments came from clustering inter-chromosomal patterns in a 4.9-billion-contact map, not from extra eigenvectors. **Symptom:** higher eigenvectors are noise, not finer biology. **Fix:** use SNIPER (inter-chr ML) or Calder (hierarchical); for differential use dcHiC.

### Hand-diffing independently phased eigenvectors
**Trigger:** subtracting/comparing two per-sample eigenvectors bin-by-bin. **Mechanism:** a weak GC correlation can flip the sign in one sample only. **Symptom:** spurious "compartment switches" concentrated on whole chromosomes/arms. **Fix:** use dcHiC (sign-coherent quantile-normalized framework) -- see hic-differential.

### Saddle strengths not comparable across samples
**Trigger:** different `n_bins`, `qrange`, resolution, or extent between compared saddles. **Mechanism:** `saddle_strength` is an extent-dependent array, not an absolute scalar. **Symptom:** strength differences that track the settings, not the biology. **Fix:** fix `n_bins`, `qrange`, resolution, and the corner extent; apply identically to all samples.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| Compartment resolution 100kb-1Mb (250kb typical) | compartment scale (Lieberman-Aiden 2009) | finer bins mix in TAD/loop structure and sparsity noise; A/B is chromosome-scale |
| Run per chromosome ARM | eigenvector-selection (cooltools docs; Mirny lab) | removes the centromere/arm gradient that otherwise hijacks E1 |
| `sort_metric='pearsonr'` | cooltools default-mismatch | default sorts by eigenvalue, so the arm gradient is reported as E1; pearsonr sorts by GC correlation |
| `n_eigs>=3` and inspect eigenvalues | eigenvector-selection | n_eigs=1 hides the arm-vs-compartment problem; which component is biology is then undeterminable |
| `clip_percentile=99.9` (cooler `eigs_cis` default) | outlier suppression | dense `cis_eig` defaults `clip_percentile=0`; the entry points differ -- do not assume |
| Saddle quantile groups ~30-50; trim 2.5% tails | cooltools tutorial | enough groups to resolve the saddle; tail trim resists outlier bins |
| Saddle strength at a fixed extent (e.g. top/bottom ~20%) | saddle_strength is an array | no canonical scalar; one extent, applied identically across samples |
| Subcompartments require deep inter-chr data | Rao 2014 (4.9B contacts) | shallow maps + extra eigenvectors give noise, not subcompartments |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `clr.matrix(balance=True)` / O/E all NaN | cooler not balanced | run `cooler balance` / `cooler.balance_cooler` first (-> matrix-operations) |
| Empty / wrong-resolution result on `.mcool` | bare `.mcool` passed | use `file.mcool::/resolutions/<bp>` URI |
| A/B compartments inverted | eigenvector sign unphased or weak phasing | pass a GC/gene-density `phasing_track` at the cooler's binning |
| E1 monotonic, no sign flips | whole-chromosome run captured the arm gradient | run per arm (`make_chromarms`); pick by max |corr| with GC |
| `frac_gc` / empty eigenvector on some chroms | chrom naming mismatch (`chr1` vs `1`) across cooler/FASTA/centromeres | harmonize names; subset the view to `clr.chromnames` |
| `saddle` KeyError on `balanced.avg` | wrong/absent expected table | pass `cooltools.expected_cis(clr, view_df=...)` output and `contact_type='cis'` |
| `AttributeError` on a cooltools function | pre-0.7 vs 0.7+ API change | `help(cooltools.eigs_cis)`; update to the viewframe signature |

## References

- Lieberman-Aiden E, van Berkum NL, et al. 2009. Comprehensive mapping of long-range interactions reveals folding principles of the human genome. *Science* 326:289-293.
- Rao SSP, Huntley MH, et al. 2014. A 3D map of the human genome at kilobase resolution reveals principles of chromatin looping. *Cell* 159:1665-1680.
- Nora EP, Goloborodko A, et al. 2017. Targeted degradation of CTCF decouples local insulation of chromosome domains from genomic compartmentalization. *Cell* 169:930-944.
- Schwarzer W, Abdennur N, et al. 2017. Two independent modes of chromatin organization revealed by cohesin removal. *Nature* 551:51-56.
- Rao SSP, Huang S-C, et al. 2017. Cohesin loss eliminates all loop domains. *Cell* 171:305-320.
- Nuebler J, Fudenberg G, Imakaev M, Abdennur N, Mirny LA. 2018. Chromatin organization by an interplay of loop extrusion and compartmental segregation. *PNAS* 115:E6697-E6706.
- Xiong K, Ma J. 2019. Revealing Hi-C subcompartments by imputing inter-chromosomal chromatin interactions. *Nat Commun* 10:5069.
- Tan L, Xing D, Chang C-H, Li H, Xie XS. 2018. Three-dimensional genome structures of single diploid human cells. *Science* 361(6405):924-928.
- Liu Y, Nanni L, et al. 2021. Systematic inference and comparison of multi-scale chromatin sub-compartments connects spatial organization to cell phenotypes. *Nat Commun* 12:2439.
- Chakraborty A, Wang JG, Ay F. 2022. dcHiC detects differential compartments across multiple Hi-C datasets. *Nat Commun* 13:6827.
- Chen Y, Zhang Y, et al. 2018. Mapping 3D genome organization relative to nuclear compartments using TSA-Seq as a cytological ruler. *J Cell Biol* 217:4025-4048.
- Abdennur N, et al. (Open2C). 2024. Cooltools: enabling high-resolution Hi-C analysis in Python. *PLoS Comput Biol* 20:e1012067.
- Abdennur N, Mirny LA. 2020. Cooler: scalable storage for Hi-C data and other genomically labeled arrays. *Bioinformatics* 36:311-316.

## Related Skills

- matrix-operations - Balancing and distance-normalized expected that compartment calling depends on
- hic-data-io - Load and access the cooler files this skill operates on
- hic-differential - dcHiC differential compartments and cross-condition switching
- tad-detection - The loop-extrusion partner of the two-mechanism framework; antagonistic strength
- hic-visualization - Render the eigenvector track and saddle plot
- chip-seq/chromatin-state-segmentation - Overlay ChromHMM/histone states on A/B compartments
- chip-seq/peak-annotation - Annotate switched bins with TF/histone peaks
- genome-intervals/bigwig-tracks - Export the eigenvector as a bigWig track
- single-cell/scatac-analysis - Single-cell chromatin context for scHi-C compartment work

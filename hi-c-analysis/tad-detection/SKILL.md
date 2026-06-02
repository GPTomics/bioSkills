---
name: bio-hi-c-analysis-tad-detection
description: Detects TAD boundaries from balanced Hi-C contact matrices via the diamond-window insulation score (cooltools insulation) and HiCExplorer hicFindTADs, returning a continuous log2 insulation track, valley-prominence boundary_strength, and Li/Otsu-thresholded is_boundary flags across a list of window sizes. Covers the multi-scale window sweep (sub-TAD to compartment-domain), why the boundary is reproducible but the domain partition is not, cross-condition comparison via differential SCORE not differential partition, and the insulation-vs-compartment orthogonality. Use when calling TADs or domain boundaries, computing insulation scores, choosing a window size, ranking boundary strength, comparing boundaries across conditions, or annotating CTCF-backed boundaries; route domain rendering to hic-visualization and boundary-feature overlap to genome-intervals.
tool_type: mixed
primary_tool: cooltools
---

## Version Compatibility

Reference examples tested with: cooler 0.10+, cooltools 0.7+, bioframe 0.7+, HiCExplorer 3.7+

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures
- CLI: `<tool> --version` then `<tool> --help` to confirm flags

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

cooltools changed its API around 0.5 -> 0.7+ (functions standardized on `view_df`/viewframe arguments; `insulation` returns the `boundary_strength_{W}`/`is_boundary_{W}` columns). A `.cool` MUST be balanced (a stored `weight` column) before insulation; `clr.matrix(balance=True)` on an unbalanced cooler returns all-NaN. A `.mcool` is multi-resolution: pass a single-resolution URI (`file.mcool::/resolutions/10000`), never the bare `.mcool`.

# TAD Detection

**"Where are the reproducible domain boundaries in my Hi-C matrix, and how strong?"** -> Compute the diamond-window insulation score on the balanced matrix, take valley minima as boundaries and their prominence as strength, and report across a LIST of window sizes rather than a single magic scale.
- Python: `cooltools.insulation(clr, [3*res, 5*res, 10*res, 25*res])` then rank by `boundary_strength_{W}`
- CLI: `hicFindTADs -m corrected.cool --outPrefix tads --correctForMultipleTesting fdr` (sweep `--minDepth/--maxDepth/--step`)

## The Single Most Important Modern Insight -- The Boundary Is Real; the Domain Is Mostly an Averaging Artifact

A population Hi-C "TAD" is the ensemble average over a heterogeneous mixture of cell-specific, stochastic domains. Single-cell imaging (Bintu 2018 *Science* 362:eaau1783) shows individual cells DO have sharp domains, but the boundary POSITION varies cell to cell - the population boundary is a *preferred* position, not a wall. Cohesin depletion abolishes population TADs while leaving single-cell domains intact, removing only the preferred-position bias. Three consequences govern every decision in this skill:

1. **"How many TADs are there" is the wrong question.** TAD number and size vary 2-5x across caller, resolution, and normalization, with NO ground truth (Forcato 2017 *Nat Methods* 14:679; Zufferey 2018 *Genome Biol* 19:217). Insulation valleys and directionality-index sign-changes give DIFFERENT boundary sets on the same matrix. The quoted "average TAD size ~880kb" is an artifact of one caller at one resolution - Zufferey explicitly states there is no average TAD size.
2. **The boundary is the reproducible, mechanistic unit.** Strong, CTCF-backed boundaries survive caller and resolution swaps; weak boundaries and exact domain extents do not. Prefer the continuous insulation/boundary-strength track over a hard domain partition for any downstream claim.
3. **The diamond window IS the analysis.** Small window (~3x bin) -> sub-TAD/fine boundaries; large window (~25x bin) -> compartment-scale domains - same matrix, totally different "TADs." Always run a list of windows and report multi-scale; never report a single-window partition as ground truth.

## TAD-Caller Taxonomy

| Method | Role | Mechanism | When |
|--------|------|-----------|------|
| cooltools `insulation` | boundary score + strength | diamond-window valleys; prominence = strength; Li/Otsu threshold | cooler-native pipeline, multi-scale, the modern default |
| HiCExplorer `hicFindTADs` | domains + boundaries + FDR | multi-window TAD-separation score with per-bin multiple-testing | CLI workflow, hierarchical sweep, FDR-controlled boundaries |
| directionality index (DI) | boundary direction | HMM on up/downstream interaction bias (Dixon 2012) | classic comparison, legacy reproducibility, gives a partition |
| Arrowhead (Juicer) | corner-score domains | arrowhead transform on the `.hic`; loop-anchored contact domains | Juicer/.hic ecosystems; calls fewer, sharper domains (Rao 2014, median ~185kb) |
| OnTAD / TADtree / rGMAP | nested/hierarchical | explicitly models meta-TADs > TADs > sub-TADs | when the question is about hierarchy or sub-TAD insulation (An 2019) |
| Stripenn / JOnTADS | stripe-aware | calls asymmetric stripes as first-class objects | when the map shows flames/stripes a flat caller mis-segments |

Insulation/DI/hicFindTADs are blind to stripes (asymmetric one-sided extrusion). Arrowhead "contact domains" are corner-anchored and categorically different from track-based boundaries - do NOT cross-compare their counts naively.

## Decision Tree by Scenario

| Scenario | Recommended | Why |
|----------|-------------|-----|
| Matrix not yet balanced | `cooler balance` / `cooler.balance_cooler` first | unbalanced insulation = coverage-driven garbage valleys |
| Reproducible boundaries, one sample | `insulation` multi-window, rank by `boundary_strength_{W}` | strength is continuous and comparable; partition is brittle |
| "What scale of domain?" | run windows `[3,5,10,25]x` bin; ~10x is the mammalian sweet spot | the window sets sub-TAD vs TAD vs compartment-domain |
| Need FDR-controlled domains | `hicFindTADs` with `--correctForMultipleTesting fdr`, sweep depths | per-bin multiple-testing on the TAD-separation score |
| Hierarchical/nested structure | OnTAD/TADtree (or compare windows) | a flat caller picks ONE level set by its window |
| Asymmetric stripes/flames present | Stripenn/JOnTADS | insulation-only pipelines are blind to stripes |
| Two conditions, boundary change | differential SCORE at matched bins, NOT intersected domain BEDs | partitions are unstable; set-differencing manufactures spurious gain/loss |
| Annotate boundaries with CTCF | -> chip-seq/peak-annotation, genome-intervals/overlap-significance | ~76-85% of boundaries are convergent CTCF + cohesin |
| Overlap boundaries with features | -> genome-intervals/interval-arithmetic | boundary BED set operations live there |
| Render domains on the matrix | -> hic-visualization | the TAD square is a colormap/resolution choice as much as a measurement |

## Insulation Score and Boundaries (multi-scale)

**Goal:** Produce a continuous boundary-strength track and threshold-flagged boundaries at several scales, so the analysis reports where insulation reproducibly dips rather than a single brittle partition.

**Approach:** Run `cooltools.insulation` on the balanced cooler with a LIST of window sizes (3-25x the bin). Each window appends its own `log2_insulation_score_{W}` (valleys = boundaries), `boundary_strength_{W}` (valley prominence - the quantitative, comparable strength), and `is_boundary_{W}` (the prominence passed through a Li histogram threshold). Rank and compare on `boundary_strength`, not on the boolean flag.

```python
import cooler
import cooltools

clr = cooler.Cooler('matrix.mcool::/resolutions/10000')   # single-resolution URI, must be balanced
res = clr.binsize
windows = [3 * res, 5 * res, 10 * res, 25 * res]   # 30k,50k,100k,250k: sub-TAD -> compartment-domain
ins = cooltools.insulation(clr, windows, verbose=True)   # clr_weight_name='weight' default -> needs ICE balancing

strong = ins[ins[f'is_boundary_{10 * res}']]   # 100kb window: ~10x bin, mammalian interphase sweet spot
ranked = ins.dropna(subset=[f'boundary_strength_{10 * res}']).sort_values(f'boundary_strength_{10 * res}', ascending=False)
```

`boundary_strength_{W}` is the scipy-style PROMINENCE of the insulation valley - continuous, quantitative, and comparable across samples; use it for ranking and cross-condition deltas. `is_boundary_{W}` is just that prominence passed through `threshold='Li'` (skimage threshold_li, an Otsu-like histogram split that is MORE PERMISSIVE than Otsu). Because the Li cutoff is fit per dataset, `is_boundary` is dataset-dependent and NOT directly comparable across samples - compare `boundary_strength`, then threshold consistently. `min_frac_valid_pixels` (default 0.66) and `min_dist_bad_bin` gate which bins get a score; sparse/blacklisted regions silently drop boundaries, so inspect `n_valid_pixels_{W}` before trusting a boundary in a low-coverage locus.

## Domains and FDR with hicFindTADs (CLI)

**Goal:** Get an FDR-controlled boundary/domain set from a multi-window TAD-separation score when a CLI workflow or hierarchical depth sweep is preferred.

**Approach:** Feed a CORRECTED (balanced) matrix and sweep the diamond depths (`--minDepth/--maxDepth/--step`); hicFindTADs computes a TAD-separation score at each depth and applies per-bin multiple-testing. The docs explicitly warn to sweep parameters before claiming a TAD count or comparing conditions.

```bash
hicFindTADs -m corrected.cool --outPrefix tads \
    --minDepth 30000 --maxDepth 100000 --step 10000 \
    --correctForMultipleTesting fdr --thresholdComparisons 0.01 --delta 0.01
# minDepth >= ~3x bin, maxDepth <= ~10x range, step >= ~2x bin; --minBoundaryDistance defaults to 4x bin
# outputs: tads_boundaries.bed, tads_domains.bed, tads_score.bedgraph, tads_tad_separation.bm, tads_zscore_matrix.h5
```

## Cross-Condition: Differential SCORE, Never Differential Partition

**Goal:** Decide which boundaries strengthen or weaken between conditions without the spurious gain/loss that comes from intersecting unstable domain calls.

**Approach:** Because partitions are unstable (caller/resolution-dependent), do NOT call TADs in each condition and set-difference the domain BEDs. Instead match resolution AND down-sample to matched valid-pixel depth, compute the bin-matched continuous insulation track at a fixed window, take the per-bin delta of `log2_insulation_score` (or `boundary_strength`), and test against a permutation/replicate null. Report boundary STRENGTHENING/WEAKENING, treating a binary boundary gain/loss as real only when strength crosses threshold robustly across replicates.

```python
ins_wt = cooltools.insulation(clr_wt, [10 * res])
ins_ko = cooltools.insulation(clr_ko, [10 * res])
key = f'log2_insulation_score_{10 * res}'
merged = ins_wt[['chrom', 'start', 'end', key]].merge(ins_ko[['chrom', 'start', 'end', key]], on=['chrom', 'start', 'end'], suffixes=('_wt', '_ko'))
merged['delta'] = merged[f'{key}_ko'] - merged[f'{key}_wt']   # negative = stronger insulation in KO; test vs a permutation null
```

Insulation (loop-extrusion barriers) and A/B compartmentalization (affinity/phase separation) are ORTHOGONAL mechanisms: CTCF degron erases insulation while compartments persist (Nora 2017 *Cell* 169:930); cohesin/RAD21 degron erases TADs+loops while compartments sharpen. Never read a boundary change as a compartment switch - a boundary can sit mid-compartment.

## Per-Method Failure Modes

### Insulation on an unbalanced matrix
**Trigger:** `insulation` on a cooler with no stored `weight` (or `clr_weight_name=None`). **Mechanism:** the diamond sum is dominated by per-bin coverage bias, not topology. **Symptom:** valleys track sequencing depth/blacklist, not domains. **Fix:** `cooler balance` first; keep the default `clr_weight_name='weight'`.

### Single magic window reported as "the TADs"
**Trigger:** calling `insulation` with one `window_bp` and treating its partition as ground truth. **Mechanism:** the window IS the scale dial; one window picks one level of a nested hierarchy. **Symptom:** sub-TAD or compartment-domain structure invisible; "TAD count" irreproducible. **Fix:** sweep `[3,5,10,25]x` bin and report multi-scale; pick the scale that matches the biological question.

### Differential partition (intersecting domain BEDs)
**Trigger:** calling TADs per condition and set-differencing the domain files. **Mechanism:** partitions are unstable, so set differences manufacture changes that are caller noise. **Symptom:** large "gained/lost TAD" lists that do not replicate. **Fix:** differential on the continuous bin-matched insulation/boundary-strength track with a permutation null.

### Window smaller than ~3x bin
**Trigger:** `window_bp < 3 * binsize`. **Mechanism:** the diamond spans too few pixels to average out noise. **Symptom:** dense spurious boundaries, no biological structure. **Fix:** set window >= 3x bin (10x is the mammalian sweet spot).

### Comparing is_boundary across samples
**Trigger:** counting `is_boundary` True in two libraries and subtracting. **Mechanism:** the Li threshold is fit per dataset; depth/strength-distribution differences shift the cutoff. **Symptom:** apparent boundary gain/loss driven by depth, not biology. **Fix:** compare continuous `boundary_strength`, then threshold consistently.

### Boundary dropped in a sparse locus
**Trigger:** a boundary expected in a low-coverage/blacklisted region is missing. **Mechanism:** `min_frac_valid_pixels` (0.66) and `min_dist_bad_bin` gate scoring; sparse diamonds get NaN. **Symptom:** no boundary where the biology predicts one. **Fix:** inspect `n_valid_pixels_{W}`; raise `min_dist_bad_bin` near bad bins or interpret cautiously.

### chrom-name mismatch (chr1 vs 1)
**Trigger:** cooler uses `chr1`, a phasing/annotation track uses `1`. **Mechanism:** chromosomes never match. **Symptom:** empty/zero output, no error. **Fix:** harmonize names across cooler, fasta, and CTCF/feature tracks.

## Quantitative Thresholds

| Threshold | Source | Rationale |
|-----------|--------|-----------|
| TAD/insulation resolution 10-40kb | domain scale | sub-Mb domains; bins must resolve boundaries without burning depth |
| Window 3-25x bin (sweep) | Open2C insulation notebook | <3x = noise; 25x = compartment-domain scale; the window is the scale dial |
| ~10x bin single window | mammalian interphase convention | e.g. 100kb window at 10kb bins for interphase TAD boundaries |
| `min_frac_valid_pixels` 0.66 | cooltools default | min valid-pixel fraction in a diamond for the bin to score |
| `threshold='Li'` | cooltools default | permissive (vs Otsu) histogram split; dataset-dependent, NOT cross-sample comparable |
| hicFindTADs depths: minDepth >=3x bin, step >=2x bin | HiCExplorer docs | the diamond depths must straddle real domain sizes; sweep before comparing |
| ~76-85% boundaries are CTCF (convergent) | Rao 2014; Vietri Rudan 2015 | strength scales with CTCF+cohesin occupancy; a sanity anchor, not a filter |
| match resolution + valid-pixel depth before gain/loss | resolution-confound | unequal depth shifts boundaries and merges sub-TADs; a false-positive engine |

## Common Errors

| Error / symptom | Cause | Solution |
|-----------------|-------|----------|
| `insulation` output all NaN | cooler not balanced | `cooler balance` / `cooler.balance_cooler` first |
| KeyError on `.mcool` / wrong resolution | bare `.mcool` passed | use `file.mcool::/resolutions/<bp>` URI |
| Dense spurious boundaries | window < ~3x bin | raise `window_bp` to >= 3x bin (10x typical) |
| Boundary counts differ wildly between samples | comparing `is_boundary` (per-dataset Li threshold) | compare continuous `boundary_strength`, threshold consistently |
| Spurious "gained/lost TADs" | differential on intersected domain partitions | differential on the continuous bin-matched score with a null |
| Empty result / missing boundary | chrom naming mismatch or sparse locus | harmonize names; inspect `n_valid_pixels_{W}` |
| `AttributeError` on cooltools call | pre-0.7 vs 0.7+ API change | `help(cooltools.insulation)`; update to the viewframe signature |

## References

- Dixon JR et al. 2012. Topological domains in mammalian genomes identified by analysis of chromatin interactions. *Nature* 485:376-380.
- Nora EP et al. 2012. Spatial partitioning of the regulatory landscape of the X-inactivation centre. *Nature* 485:381-385.
- Crane E et al. 2015. Condensin-driven remodelling of X chromosome topology during dosage compensation. *Nature* 523:240-244. (Introduced the diamond-window insulation score.)
- Rao SSP et al. 2014. A 3D map of the human genome at kilobase resolution reveals principles of chromatin looping. *Cell* 159:1665-1680. (Arrowhead contact domains; convergent CTCF.)
- Fudenberg G et al. 2016. Formation of chromosomal domains by loop extrusion. *Cell Rep* 15:2038-2049.
- Forcato M et al. 2017. Comparison of computational methods for Hi-C data analysis. *Nat Methods* 14:679-685.
- Nora EP et al. 2017. Targeted degradation of CTCF decouples local insulation of chromosome domains from genomic compartmentalization. *Cell* 169:930-944.
- Bintu B et al. 2018. Super-resolution chromatin tracing reveals domains and cooperative interactions in single cells. *Science* 362:eaau1783.
- Vian L et al. 2018. The energetics and physiological impact of cohesin extrusion. *Cell* 173:1165-1178. (Architectural stripes from one-sided extrusion.)
- Zufferey M, Tavernari D, Oricchio E, Ciriello G. 2018. Comparison of computational methods for the identification of topologically associating domains. *Genome Biol* 19:217.
- An L, Yang T, Yang J et al. 2019. OnTAD: hierarchical domain structure reveals the divergence of activity among TADs and boundaries. *Genome Biol* 20:282.
- Lupianez DG et al. 2015. Disruptions of topological chromatin domains cause pathogenic rewiring of gene-enhancer interactions. *Cell* 161:1012-1025.
- Vietri Rudan M, Barrington C, Henderson S et al. 2015. Comparative Hi-C reveals that CTCF underlies evolution of chromosomal domain architecture. *Cell Rep* 10(8):1297-1309.
- Open2C, Abdennur N, Abraham S, Fudenberg G, Flyamer IM, Galitsyna AA et al. 2024. Cooltools: enabling high-resolution Hi-C analysis in Python. *PLoS Comput Biol* 20:e1012067.
- Ramirez F et al. 2018. High-resolution TADs reveal DNA sequences underlying genome organization in cells. *Nat Commun* 9:189. (HiCExplorer.)

## Related Skills

- matrix-operations - Balancing and O/E that insulation scoring depends on
- hic-data-io - Load and access the cooler files this skill operates on
- compartment-analysis - The orthogonal Mb-scale mechanism; a boundary is not a compartment switch
- loop-calling - Convergent-CTCF loops anchor the strongest boundaries; stripes need a stripe-aware caller
- hic-differential - Replicate-aware cross-condition contact comparison
- hic-visualization - Render domains/boundaries on the contact matrix
- chip-seq/peak-annotation - Annotate boundaries with CTCF/cohesin peaks
- genome-intervals/interval-arithmetic - Overlap boundary BEDs with features
- genome-intervals/overlap-significance - Test boundary/CTCF co-localization against a matched null
